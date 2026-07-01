#!/usr/bin/env node

import fs from "node:fs";
import path from "node:path";
import readline from "node:readline";

const cwd = process.cwd();
const resultsDir = path.join(cwd, "results");
const configPath = path.join(cwd, "examples", "trp_l_tryptophan.json");
const dumpPath = path.join(resultsDir, "trp_l_tryptophan.dump");
const logPath = path.join(resultsDir, "trp_l_tryptophan.log");
const ovitoSurfacePath = path.join(resultsDir, "trp_l_tryptophan_ovito_surface_area.csv");
const surfaceEnergyPath = path.join(resultsDir, "trp_l_tryptophan_surface_energy_ovito.csv");
const finalPointsPath = path.join(resultsDir, "trp_l_tryptophan_final_points.tsv");
const finalSurfacePointsPath = path.join(resultsDir, "trp_l_tryptophan_final_surface_points.tsv");

const evPerA2ToMJPerM2 = 16021.76634;

function assertFile(filePath) {
  if (!fs.existsSync(filePath)) {
    throw new Error(`Missing required file: ${filePath}`);
  }
}

function parseLog() {
  assertFile(logPath);
  return fs.readFileSync(logPath, "utf8").trim().split(/\r?\n/).map((line, index) => {
    const time = Number(line.match(/TIME=([^\s]+)/)?.[1]);
    const energy = Number(line.match(/ENERGY=([^\s]+)/)?.[1]);
    const occupancy = Number(line.match(/OCCUPANCY=([^\s]+)/)?.[1]);
    if (!Number.isFinite(time) || !Number.isFinite(energy) || !Number.isFinite(occupancy)) {
      throw new Error(`Could not parse log line ${index + 1}: ${line}`);
    }
    return { time, energy, occupancy };
  });
}

function parseCsv(filePath) {
  assertFile(filePath);
  const [headerLine, ...lines] = fs.readFileSync(filePath, "utf8").trim().split(/\r?\n/);
  const headers = headerLine.split(",");
  return lines.map((line) => {
    const values = line.split(",");
    return Object.fromEntries(headers.map((header, index) => [header, values[index]]));
  });
}

function loadConfig() {
  assertFile(configPath);
  return JSON.parse(fs.readFileSync(configPath, "utf8"));
}

function positiveModulo(value, modulus) {
  return ((value % modulus) + modulus) % modulus;
}

function interactionEnergy(distance, interactions) {
  let index = 0;
  while (distance > interactions.cutoffs[index]) {
    index += 1;
  }
  return interactions.interaction_energies[index];
}

function neighborTables(config) {
  const lattice = config.lattice;
  const interactions = config.interactions;
  const [a, b, c] = lattice.lattice_parameters;
  const basis = lattice.atomic_basis;
  const maxCutoff = interactions.cutoffs.at(-1);
  const maxCellOffset = [
    Math.ceil(maxCutoff / a) + 1,
    Math.ceil(maxCutoff / b) + 1,
    Math.ceil(maxCutoff / c) + 1,
  ];
  const tables = Array.from({ length: basis.length }, () => []);
  for (let fromBasis = 0; fromBasis < basis.length; fromBasis += 1) {
    for (let toBasis = 0; toBasis < basis.length; toBasis += 1) {
      for (let di = -maxCellOffset[0]; di <= maxCellOffset[0]; di += 1) {
        for (let dj = -maxCellOffset[1]; dj <= maxCellOffset[1]; dj += 1) {
          for (let dk = -maxCellOffset[2]; dk <= maxCellOffset[2]; dk += 1) {
            const dx = (di + basis[toBasis][0] - basis[fromBasis][0]) * a;
            const dy = (dj + basis[toBasis][1] - basis[fromBasis][1]) * b;
            const dz = (dk + basis[toBasis][2] - basis[fromBasis][2]) * c;
            const distance = Math.hypot(dx, dy, dz);
            if (distance > 0.0 && distance <= maxCutoff) {
              tables[fromBasis].push({
                di,
                dj,
                dk,
                toBasis,
                energy: interactionEnergy(distance, interactions),
              });
            }
          }
        }
      }
    }
  }
  return tables;
}

function siteHelpers(config) {
  const [nx, ny, nz] = config.lattice.dimensions;
  const basisCount = config.lattice.atomic_basis.length;
  return {
    idToCell(id) {
      const siteBasis = id % basisCount;
      const cell = Math.floor(id / basisCount);
      const k = cell % nz;
      const j = Math.floor(cell / nz) % ny;
      const i = Math.floor(cell / (ny * nz));
      return { i, j, k, siteBasis };
    },
    cellToId(i, j, k, siteBasis) {
      return (((i * ny + j) * nz + k) * basisCount) + siteBasis;
    },
  };
}

function surfaceInfo(ids, occupied, neighbors, config, wantMask = false) {
  const [nx, ny, nz] = config.lattice.dimensions;
  const { idToCell, cellToId } = siteHelpers(config);
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

function cohesiveEnergy(config, neighbors) {
  const totalPatternEnergy = neighbors.flat().reduce((sum, neighbor) => sum + neighbor.energy, 0.0);
  return 0.5 * totalPatternEnergy / config.lattice.atomic_basis.length;
}

function density(config) {
  const [a, b, c] = config.lattice.lattice_parameters;
  return config.lattice.atomic_basis.length / (a * b * c);
}

function initialSurfaceDensity(config, neighbors) {
  const [nx, ny, nz] = config.lattice.dimensions;
  const [a, b, c] = config.lattice.lattice_parameters;
  const basis = config.lattice.atomic_basis;
  const basisCount = basis.length;
  const siteCount = nx * ny * nz * basisCount;
  const centerBasis = basis.reduce((sum, item) => [
    sum[0] + item[0],
    sum[1] + item[1],
    sum[2] + item[2],
  ], [0, 0, 0]).map((value) => value / basisCount);
  const center = [
    ((nx - 1) * 0.5 + centerBasis[0]) * a,
    ((ny - 1) * 0.5 + centerBasis[1]) * b,
    ((nz - 1) * 0.5 + centerBasis[2]) * c,
  ];
  const occupied = new Set();
  const radiusSquared = config.growth.initial_radius ** 2;
  const { idToCell } = siteHelpers(config);
  for (let id = 0; id < siteCount; id += 1) {
    const { i, j, k, siteBasis } = idToCell(id);
    const position = [
      (i + basis[siteBasis][0]) * a,
      (j + basis[siteBasis][1]) * b,
      (k + basis[siteBasis][2]) * c,
    ];
    const dx = position[0] - center[0];
    const dy = position[1] - center[1];
    const dz = position[2] - center[2];
    if (dx * dx + dy * dy + dz * dz <= radiusSquared) {
      occupied.add(id);
    }
  }
  const ids = Int32Array.from(occupied);
  const solidSurfaceSites = surfaceInfo(ids, occupied, neighbors, config).count;
  const initialSurfaceArea = 4.0 * Math.PI * config.growth.initial_radius ** 2;
  return solidSurfaceSites / initialSurfaceArea;
}

function kappa(config, neighbors) {
  const rho = density(config);
  const surfaceDensity = initialSurfaceDensity(config, neighbors);
  return ((4.0 / 3.0 * Math.PI * rho) ** (1.0 / 3.0)
      * config.solvent.diffusivity
      * config.solvent.solubility_limit)
    / (surfaceDensity * config.growth.desired_size ** (1.0 / 3.0));
}

async function parseDump(logEntries, neighbors, config) {
  assertFile(dumpPath);
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
    const info = surfaceInfo(ids, occupied, neighbors, config, wantMask);
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
      };
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

function writeSurfaceEnergyCsv(logEntries, surfaceRows, config, neighbors, frames) {
  const totalSites = config.lattice.dimensions.reduce((product, value) => product * value, 1)
    * config.lattice.atomic_basis.length;
  const ecoh = cohesiveEnergy(config, neighbors);
  const timeScale = kappa(config, neighbors);

  if (surfaceRows.length !== logEntries.length) {
    throw new Error(`Frame/log mismatch: OVITO CSV has ${surfaceRows.length} rows, log has ${logEntries.length}`);
  }

  const rows = [[
    "frame",
    "step",
    "time_s",
    "time_kappa",
    "energy_eV",
    "n_atoms",
    "ovito_surface_area_A2",
    "gamma_ovito_eV_A2",
    "gamma_ovito_mJ_m2",
  ].join(",")];

  for (const row of surfaceRows) {
    const frame = Number(row.frame);
    const log = logEntries[frame];
    const surfaceArea = Number(row.surface_area_A2);
    const atomsFromOvito = Number(row.n_atoms);
    const nAtoms = Number.isFinite(atomsFromOvito)
      ? atomsFromOvito
      : Math.round(log.occupancy * totalSites);
    const step = row.timestep && row.timestep !== "undefined" ? row.timestep : frames[frame].step;
    const gammaEVPerA2 = (log.energy - ecoh * nAtoms) / surfaceArea;
    rows.push([
      frame,
      step,
      log.time,
      log.time * timeScale,
      log.energy,
      nAtoms,
      surfaceArea,
      gammaEVPerA2,
      gammaEVPerA2 * evPerA2ToMJPerM2,
    ].join(","));
  }

  fs.writeFileSync(surfaceEnergyPath, `${rows.join("\n")}\n`);
  console.error(`Wrote ${surfaceEnergyPath}`);
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

  fs.writeFileSync(finalPointsPath, `${allRows.join("\n")}\n`);
  fs.writeFileSync(finalSurfacePointsPath, `${surfaceRows.join("\n")}\n`);
  console.error(`Wrote ${finalPointsPath}`);
  console.error(`Wrote ${finalSurfacePointsPath}`);
}

async function main() {
  const config = loadConfig();
  const neighbors = neighborTables(config);
  const logEntries = parseLog();
  const surfaceRows = parseCsv(ovitoSurfacePath);
  const { frames, finalFrame } = await parseDump(logEntries, neighbors, config);

  writeSurfaceEnergyCsv(logEntries, surfaceRows, config, neighbors, frames);
  writeFinalPointClouds(finalFrame);
  console.error(`Cohesive energy: ${cohesiveEnergy(config, neighbors).toFixed(6)} eV/molecule`);
  console.error(`kappa: ${kappa(config, neighbors).toFixed(6)} 1/s`);
}

main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
