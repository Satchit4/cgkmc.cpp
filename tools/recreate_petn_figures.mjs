#!/usr/bin/env node

import fs from "node:fs";
import path from "node:path";
import readline from "node:readline";

const repoRoot = process.cwd();
const resultsDir = path.join(repoRoot, "results");
const dumpPath = path.join(resultsDir, "petn.dump");
const logPath = path.join(resultsDir, "petn.log");

const dims = [62, 62, 140];
const latticeParameters = [9.087, 9.087, 6.738];
const basis = [
  [0.0, 0.0, 0.0],
  [0.5, 0.5, 0.5],
];
const cutoffs = [7.0, 7.5, 9.5];
const initialRadius = 75.0;
const desiredSize = 40000;
const diffusivity = 1.0e11;
const solubilityLimit = 1.0e-4;
const cohesiveEnergyEV = -1.034;
const aeSurfaceEnergyMJm2 = 88.0;
const evPerA2ToMJPerM2 = 16021.76634;

const [nx, ny, nz] = dims;
const [a, b, c] = latticeParameters;
const basisCount = basis.length;
const totalSites = nx * ny * nz * basisCount;
const density = basisCount / (a * b * c);

function assertFile(filePath) {
  if (!fs.existsSync(filePath)) {
    throw new Error(`Missing required file: ${filePath}`);
  }
}

function parseLog() {
  const text = fs.readFileSync(logPath, "utf8").trim();
  if (!text) {
    throw new Error(`Empty log file: ${logPath}`);
  }
  return text.split(/\r?\n/).map((line, index) => {
    const time = Number(line.match(/TIME=([^\s]+)/)?.[1]);
    const energy = Number(line.match(/ENERGY=([^\s]+)/)?.[1]);
    const occupancy = Number(line.match(/OCCUPANCY=([^\s]+)/)?.[1]);
    if (!Number.isFinite(time) || !Number.isFinite(energy) || !Number.isFinite(occupancy)) {
      throw new Error(`Could not parse log line ${index + 1}: ${line}`);
    }
    return { time, energy, occupancy };
  });
}

function neighborTables() {
  const maxCutoff = cutoffs.at(-1);
  const maxCellOffset = [
    Math.ceil(maxCutoff / a) + 1,
    Math.ceil(maxCutoff / b) + 1,
    Math.ceil(maxCutoff / c) + 1,
  ];
  const tables = [[], []];
  for (let fromBasis = 0; fromBasis < basisCount; fromBasis += 1) {
    for (let toBasis = 0; toBasis < basisCount; toBasis += 1) {
      for (let di = -maxCellOffset[0]; di <= maxCellOffset[0]; di += 1) {
        for (let dj = -maxCellOffset[1]; dj <= maxCellOffset[1]; dj += 1) {
          for (let dk = -maxCellOffset[2]; dk <= maxCellOffset[2]; dk += 1) {
            const dx = (di + basis[toBasis][0] - basis[fromBasis][0]) * a;
            const dy = (dj + basis[toBasis][1] - basis[fromBasis][1]) * b;
            const dz = (dk + basis[toBasis][2] - basis[fromBasis][2]) * c;
            const dist = Math.hypot(dx, dy, dz);
            if (dist > 1.0e-9 && dist <= maxCutoff + 1.0e-9) {
              tables[fromBasis].push({ di, dj, dk, toBasis });
            }
          }
        }
      }
    }
  }
  return tables;
}

function positiveModulo(value, modulus) {
  return ((value % modulus) + modulus) % modulus;
}

function idToCell(id) {
  const siteBasis = id % basisCount;
  const cell = Math.floor(id / basisCount);
  const k = cell % nz;
  const j = Math.floor(cell / nz) % ny;
  const i = Math.floor(cell / (ny * nz));
  return { i, j, k, siteBasis };
}

function cellToId(i, j, k, siteBasis) {
  return (((i * ny + j) * nz + k) * basisCount) + siteBasis;
}

function surfaceInfo(ids, occupied, neighbors, wantMask = false) {
  let count = 0;
  const mask = wantMask ? new Uint8Array(ids.length) : null;
  for (let index = 0; index < ids.length; index += 1) {
    const id = ids[index];
    const { i, j, k, siteBasis } = idToCell(id);
    let exposed = false;
    for (const neighbor of neighbors[siteBasis]) {
      const ni = positiveModulo(i + neighbor.di, nx);
      const nj = positiveModulo(j + neighbor.dj, ny);
      const nk = positiveModulo(k + neighbor.dk, nz);
      const neighborId = cellToId(ni, nj, nk, neighbor.toBasis);
      if (!occupied.has(neighborId)) {
        exposed = true;
        break;
      }
    }
    if (exposed) {
      count += 1;
      if (mask) {
        mask[index] = 1;
      }
    }
  }
  return { count, mask };
}

async function parseDump(logEntries, neighbors) {
  const stream = fs.createReadStream(dumpPath);
  const rl = readline.createInterface({ input: stream, crlfDelay: Infinity });

  const frames = [];
  let frame = null;
  let state = "seek";
  let skipBounds = 0;
  let remainingAtoms = 0;
  let atomIndex = 0;
  let ids = null;
  let xs = null;
  let ys = null;
  let zs = null;
  let finalFrame = null;

  function processFrame() {
    const occupied = new Set(ids);
    const wantMask = frames.length === logEntries.length - 1;
    const info = surfaceInfo(ids, occupied, neighbors, wantMask);
    frames.push({
      step: frame.step,
      dumpTime: frame.dumpTime,
      atoms: frame.atoms,
      surfaceSites: info.count,
    });
    if (wantMask) {
      finalFrame = {
        ids,
        xs,
        ys,
        zs,
        surfaceMask: info.mask,
        step: frame.step,
      };
    }
    if (frames.length % 100 === 0 || frames.length === 1) {
      console.error(`Processed ${frames.length} / ${logEntries.length} frames`);
    }
  }

  for await (const rawLine of rl) {
    const line = rawLine.trim();
    if (line === "ITEM: TIMESTEP") {
      frame = {};
      state = "timestep";
      continue;
    }
    if (state === "timestep") {
      const [stepString, timeString] = line.split(/\s+/);
      frame.step = Number(stepString);
      frame.dumpTime = Number(timeString);
      state = "numberHeader";
      continue;
    }
    if (line === "ITEM: NUMBER OF ATOMS") {
      state = "atomCount";
      continue;
    }
    if (state === "atomCount") {
      frame.atoms = Number(line);
      ids = new Int32Array(frame.atoms);
      const isFinal = frames.length === logEntries.length - 1;
      xs = isFinal ? new Float64Array(frame.atoms) : null;
      ys = isFinal ? new Float64Array(frame.atoms) : null;
      zs = isFinal ? new Float64Array(frame.atoms) : null;
      state = "boundsHeader";
      continue;
    }
    if (line === "ITEM: BOX BOUNDS xy xz xx yy zz") {
      skipBounds = 3;
      state = "bounds";
      continue;
    }
    if (state === "bounds" && skipBounds > 0) {
      skipBounds -= 1;
      continue;
    }
    if (line === "ITEM: ATOMS id type x y z") {
      remainingAtoms = frame.atoms;
      atomIndex = 0;
      state = "atoms";
      continue;
    }
    if (state === "atoms" && remainingAtoms > 0) {
      const parts = line.split(/\s+/);
      ids[atomIndex] = Number(parts[0]);
      if (xs) {
        xs[atomIndex] = Number(parts[2]);
        ys[atomIndex] = Number(parts[3]);
        zs[atomIndex] = Number(parts[4]);
      }
      atomIndex += 1;
      remainingAtoms -= 1;
      if (remainingAtoms === 0) {
        processFrame();
        state = "seek";
      }
    }
  }

  if (frames.length !== logEntries.length) {
    throw new Error(`Frame/log mismatch: dump has ${frames.length}, log has ${logEntries.length}`);
  }
  return { frames, finalFrame };
}

function writeCsv(logEntries, frames) {
  const initialSurfaceArea = 4.0 * Math.PI * initialRadius ** 2;
  const surfaceDensity = frames[0].surfaceSites / initialSurfaceArea;
  const kappa = ((4.0 / 3.0 * Math.PI * density) ** (1.0 / 3.0) * diffusivity * solubilityLimit) /
    (surfaceDensity * desiredSize ** (1.0 / 3.0));

  const rows = [
    [
      "frame",
      "step",
      "time_s",
      "time_kappa",
      "energy_eV",
      "n_atoms",
      "surface_sites",
      "surface_area_A2",
      "gamma_eV_A2",
      "gamma_mJ_m2",
      "spherical_surface_area_A2",
      "gamma_spherical_mJ_m2",
    ].join(","),
  ];
  const data = frames.map((frame, index) => {
    const log = logEntries[index];
    const surfaceArea = frame.surfaceSites / surfaceDensity;
    const gammaEVPerA2 = (log.energy - cohesiveEnergyEV * frame.atoms) / surfaceArea;
    const gammaMJPerM2 = gammaEVPerA2 * evPerA2ToMJPerM2;
    const equivalentRadius = (3.0 * frame.atoms / (4.0 * Math.PI * density)) ** (1.0 / 3.0);
    const sphericalSurfaceArea = 4.0 * Math.PI * equivalentRadius ** 2;
    const gammaSphericalMJPerM2 = ((log.energy - cohesiveEnergyEV * frame.atoms) / sphericalSurfaceArea) *
      evPerA2ToMJPerM2;
    const item = {
      frame: index,
      step: frame.step,
      time: log.time,
      timeKappa: log.time * kappa,
      energy: log.energy,
      atoms: frame.atoms,
      surfaceSites: frame.surfaceSites,
      surfaceArea,
      gammaEVPerA2,
      gammaMJPerM2,
      sphericalSurfaceArea,
      gammaSphericalMJPerM2,
    };
    rows.push([
      item.frame,
      item.step,
      item.time,
      item.timeKappa,
      item.energy,
      item.atoms,
      item.surfaceSites,
      item.surfaceArea,
      item.gammaEVPerA2,
      item.gammaMJPerM2,
      item.sphericalSurfaceArea,
      item.gammaSphericalMJPerM2,
    ].join(","));
    return item;
  });

  fs.writeFileSync(path.join(resultsDir, "petn_surface_energy.csv"), `${rows.join("\n")}\n`);
  return { data, surfaceDensity, kappa };
}

function svgEscape(text) {
  return String(text)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;");
}

function polylinePath(points, sx, sy) {
  return points.map((point, index) => {
    const command = index === 0 ? "M" : "L";
    return `${command}${sx(point.timeKappa).toFixed(2)},${sy(point.gammaMJPerM2).toFixed(2)}`;
  }).join(" ");
}

function figure2Svg(data, options = {}) {
  const width = 1280;
  const height = 900;
  const margin = { top: 70, right: 55, bottom: 100, left: 120 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  const maxX = options.maxX ?? Math.ceil(Math.max(...data.map((item) => item.timeKappa)));
  const plotData = data.filter((item) => item.timeKappa <= maxX);
  const minY = 80;
  const maxY = Math.ceil(Math.max(210, Math.max(...plotData.map((item) => item.gammaMJPerM2)) + 5) / 10) * 10;
  const sx = (x) => margin.left + (x / maxX) * plotWidth;
  const sy = (y) => margin.top + (1.0 - (y - minY) / (maxY - minY)) * plotHeight;

  const grid = [];
  for (let x = 0; x <= maxX; x += 5) {
    grid.push(`<line x1="${sx(x)}" y1="${margin.top}" x2="${sx(x)}" y2="${margin.top + plotHeight}" class="grid"/>`);
    grid.push(`<text x="${sx(x)}" y="${height - 58}" class="tick" text-anchor="middle">${x}</text>`);
  }
  for (let y = minY; y <= maxY; y += 20) {
    grid.push(`<line x1="${margin.left}" y1="${sy(y)}" x2="${margin.left + plotWidth}" y2="${sy(y)}" class="grid"/>`);
    grid.push(`<text x="${margin.left - 18}" y="${sy(y) + 8}" class="tick" text-anchor="end">${y}</text>`);
  }

  const line = polylinePath(plotData, sx, sy);
  const aeY = sy(aeSurfaceEnergyMJm2);

  return `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">
  <style>
    .axis { stroke: #111; stroke-width: 2.2; fill: none; }
    .grid { stroke: #bdbdbd; stroke-width: 1; opacity: 0.75; }
    .tick { font: 24px Arial, sans-serif; fill: #111; }
    .label { font: 30px Arial, sans-serif; fill: #111; }
    .legend { font: 23px Arial, sans-serif; fill: #111; }
  </style>
  <rect width="100%" height="100%" fill="white"/>
  ${grid.join("\n  ")}
  <path d="${line}" fill="none" stroke="#ff7f0e" stroke-width="3.2" stroke-linejoin="round" stroke-linecap="round"/>
  <line x1="${margin.left}" y1="${aeY}" x2="${margin.left + plotWidth}" y2="${aeY}" stroke="#111" stroke-width="2.6" stroke-dasharray="10 8"/>
  <rect x="${margin.left}" y="${margin.top}" width="${plotWidth}" height="${plotHeight}" class="axis"/>
  <text x="${margin.left + plotWidth / 2}" y="${height - 18}" class="label" text-anchor="middle">time (units of kappa^-1)</text>
  <text x="38" y="${margin.top + plotHeight / 2}" class="label" text-anchor="middle" transform="rotate(-90 38 ${margin.top + plotHeight / 2})">surface energy (mJ/m^2)</text>
  <rect x="${width - 430}" y="${margin.top + 18}" width="335" height="105" rx="4" fill="white" stroke="#d0d0d0"/>
  <line x1="${width - 405}" y1="${margin.top + 52}" x2="${width - 330}" y2="${margin.top + 52}" stroke="#ff7f0e" stroke-width="3.2"/>
  <text x="${width - 310}" y="${margin.top + 60}" class="legend">KMC simulation</text>
  <line x1="${width - 405}" y1="${margin.top + 95}" x2="${width - 330}" y2="${margin.top + 95}" stroke="#111" stroke-width="2.6" stroke-dasharray="10 8"/>
  <text x="${width - 310}" y="${margin.top + 103}" class="legend">AE model</text>
</svg>
`;
}

function rotatePoint(point, orientation) {
  const [x0, y0, z0] = point;
  const cy = Math.cos(orientation.yaw);
  const sy = Math.sin(orientation.yaw);
  const cp = Math.cos(orientation.pitch);
  const sp = Math.sin(orientation.pitch);
  const cr = Math.cos(orientation.roll);
  const sr = Math.sin(orientation.roll);

  const x1 = cy * x0 - sy * y0;
  const y1 = sy * x0 + cy * y0;
  const z1 = z0;

  const x2 = x1;
  const y2 = cp * y1 - sp * z1;
  const z2 = sp * y1 + cp * z1;

  return [
    cr * x2 - sr * z2,
    y2,
    sr * x2 + cr * z2,
  ];
}

function collectPoints(finalFrame, surfaceOnly) {
  const points = [];
  let cx = 0;
  let cy = 0;
  let cz = 0;
  for (let i = 0; i < finalFrame.ids.length; i += 1) {
    cx += finalFrame.xs[i];
    cy += finalFrame.ys[i];
    cz += finalFrame.zs[i];
  }
  cx /= finalFrame.ids.length;
  cy /= finalFrame.ids.length;
  cz /= finalFrame.ids.length;

  for (let i = 0; i < finalFrame.ids.length; i += 1) {
    if (surfaceOnly && !finalFrame.surfaceMask[i]) {
      continue;
    }
    points.push([
      finalFrame.xs[i] - cx,
      finalFrame.ys[i] - cy,
      finalFrame.zs[i] - cz,
    ]);
  }
  return points;
}

function writeFinalPointClouds(finalFrame) {
  let cx = 0;
  let cy = 0;
  let cz = 0;
  for (let i = 0; i < finalFrame.ids.length; i += 1) {
    cx += finalFrame.xs[i];
    cy += finalFrame.ys[i];
    cz += finalFrame.zs[i];
  }
  cx /= finalFrame.ids.length;
  cy /= finalFrame.ids.length;
  cz /= finalFrame.ids.length;

  const allRows = ["id\tx\ty\tz"];
  const surfaceRows = ["id\tx\ty\tz"];
  for (let i = 0; i < finalFrame.ids.length; i += 1) {
    const row = [
      finalFrame.ids[i],
      finalFrame.xs[i] - cx,
      finalFrame.ys[i] - cy,
      finalFrame.zs[i] - cz,
    ].join("\t");
    allRows.push(row);
    if (finalFrame.surfaceMask[i]) {
      surfaceRows.push(row);
    }
  }

  fs.writeFileSync(path.join(resultsDir, "petn_final_points.tsv"), `${allRows.join("\n")}\n`);
  fs.writeFileSync(path.join(resultsDir, "petn_final_surface_points.tsv"), `${surfaceRows.join("\n")}\n`);
}

function projected(points, orientation) {
  return points.map((point) => {
    const [x, depth, z] = rotatePoint(point, orientation);
    return { x, y: -z, depth };
  });
}

function panelGeometry(projectedPoints, panel) {
  const xs = projectedPoints.map((point) => point.x);
  const ys = projectedPoints.map((point) => point.y);
  const minX = Math.min(...xs);
  const maxX = Math.max(...xs);
  const minY = Math.min(...ys);
  const maxY = Math.max(...ys);
  const scale = Math.min(panel.width / (maxX - minX), panel.height / (maxY - minY)) * 0.9;
  const offsetX = panel.x + panel.width / 2 - ((minX + maxX) / 2) * scale;
  const offsetY = panel.y + panel.height / 2 - ((minY + maxY) / 2) * scale;
  return { scale, offsetX, offsetY };
}

function depthColor(base, depth, minDepth, maxDepth) {
  const t = (depth - minDepth) / Math.max(1.0e-9, maxDepth - minDepth);
  const shade = 0.76 + 0.24 * t;
  return `rgb(${Math.round(base[0] * shade)}, ${Math.round(base[1] * shade)}, ${Math.round(base[2] * shade)})`;
}

function drawAxisTriad(orientation, x, y, scale) {
  const axes = [
    { vector: [1, 0, 0], label: "(100)" },
    { vector: [0, 1, 0], label: "(010)" },
    { vector: [0, 0, 1], label: "(001)" },
  ];
  const lines = [];
  for (const axis of axes) {
    const [rx, , rz] = rotatePoint(axis.vector, orientation);
    const x2 = x + rx * scale;
    const y2 = y - rz * scale;
    lines.push(`<line x1="${x}" y1="${y}" x2="${x2}" y2="${y2}" stroke="#111" stroke-width="2" marker-end="url(#arrow)"/>`);
    lines.push(`<text x="${x2 + Math.sign(rx || 1) * 6}" y="${y2 + Math.sign(-rz || 1) * 13}" class="axisLabel">${axis.label}</text>`);
  }
  return lines.join("\n");
}

function pointElements(projectedPoints, panel, options) {
  const geometry = panelGeometry(projectedPoints, panel);
  const depths = projectedPoints.map((point) => point.depth);
  const minDepth = Math.min(...depths);
  const maxDepth = Math.max(...depths);
  const sorted = projectedPoints.toSorted((left, right) => left.depth - right.depth);
  const elements = [];
  const stride = options.stride ?? 1;
  for (let i = 0; i < sorted.length; i += stride) {
    const point = sorted[i];
    const cx = geometry.offsetX + point.x * geometry.scale;
    const cy = geometry.offsetY + point.y * geometry.scale;
    const fill = depthColor(options.baseColor, point.depth, minDepth, maxDepth);
    elements.push(`<circle cx="${cx.toFixed(2)}" cy="${cy.toFixed(2)}" r="${options.radius}" fill="${fill}" opacity="${options.opacity}"/>`);
  }
  return elements.join("\n");
}

function figure3Svg(finalFrame) {
  const width = 1500;
  const height = 1120;
  const panels = [
    { x: 90, y: 55, width: 590, height: 430 },
    { x: 815, y: 55, width: 590, height: 430 },
    { x: 90, y: 585, width: 590, height: 430 },
    { x: 815, y: 585, width: 590, height: 430 },
  ];
  const orientations = [
    { yaw: -0.65, pitch: -0.95, roll: 0.05 },
    { yaw: 0.58, pitch: -1.02, roll: -0.2 },
  ];
  const allPoints = collectPoints(finalFrame, false);
  const surfacePoints = collectPoints(finalFrame, true);
  const projectedAllA = projected(allPoints, orientations[0]);
  const projectedAllB = projected(allPoints, orientations[1]);
  const projectedSurfaceA = projected(surfacePoints, orientations[0]);
  const projectedSurfaceB = projected(surfacePoints, orientations[1]);

  return `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">
  <defs>
    <marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto" markerUnits="strokeWidth">
      <path d="M0,0 L8,4 L0,8 z" fill="#111"/>
    </marker>
  </defs>
  <style>
    .axisLabel { font: 15px Arial, sans-serif; fill: #111; }
    .facet { font: 26px Arial, sans-serif; fill: #161016; opacity: 0.85; }
    .caption { font: 22px Arial, sans-serif; fill: #111; }
  </style>
  <rect width="100%" height="100%" fill="white"/>
  ${pointElements(projectedAllA, panels[0], { baseColor: [230, 72, 70], radius: 0.95, opacity: 0.42, stride: 1 })}
  ${pointElements(projectedAllB, panels[1], { baseColor: [230, 72, 70], radius: 0.95, opacity: 0.42, stride: 1 })}
  ${pointElements(projectedSurfaceA, panels[2], { baseColor: [164, 73, 205], radius: 2.3, opacity: 0.78, stride: 1 })}
  ${pointElements(projectedSurfaceB, panels[3], { baseColor: [164, 73, 205], radius: 2.3, opacity: 0.78, stride: 1 })}
  ${drawAxisTriad(orientations[0], panels[0].x + 85, panels[0].y + panels[0].height - 70, 45)}
  ${drawAxisTriad(orientations[1], panels[1].x + 85, panels[1].y + panels[1].height - 70, 45)}
  ${drawAxisTriad(orientations[0], panels[2].x + 85, panels[2].y + panels[2].height - 70, 45)}
  ${drawAxisTriad(orientations[1], panels[3].x + 85, panels[3].y + panels[3].height - 70, 45)}
  <text x="${panels[2].x + 260}" y="${panels[2].y + 95}" class="facet">(101)</text>
  <text x="${panels[2].x + 335}" y="${panels[2].y + 215}" class="facet">(110)</text>
  <text x="${panels[2].x + 185}" y="${panels[2].y + 265}" class="facet">(1\u030510)</text>
  <text x="${panels[2].x + 290}" y="${panels[2].y + 350}" class="facet">(10\u03051)</text>
  <text x="${panels[3].x + 300}" y="${panels[3].y + 170}" class="facet">(1\u030510)</text>
  <text x="${panels[3].x + 185}" y="${panels[3].y + 260}" class="facet">(01\u03051)</text>
  <text x="${width / 2}" y="${height - 32}" class="caption" text-anchor="middle">Final PETN morphology from the last cgkmc dump frame, rendered in two orientations.</text>
</svg>
`;
}

async function main() {
  assertFile(dumpPath);
  assertFile(logPath);
  const logEntries = parseLog();
  const neighbors = neighborTables();
  console.error(`Neighbor counts by basis: ${neighbors.map((items) => items.length).join(", ")}`);
  console.error(`Total lattice sites: ${totalSites}`);
  const { frames, finalFrame } = await parseDump(logEntries, neighbors);
  writeFinalPointClouds(finalFrame);

  fs.writeFileSync(path.join(resultsDir, "figure3_petn_final_morphology.svg"), figure3Svg(finalFrame));

  console.error(`Wrote ${path.join(resultsDir, "petn_final_points.tsv")}`);
  console.error(`Wrote ${path.join(resultsDir, "petn_final_surface_points.tsv")}`);
  console.error(`Wrote ${path.join(resultsDir, "figure3_petn_final_morphology.svg")}`);
}

main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
