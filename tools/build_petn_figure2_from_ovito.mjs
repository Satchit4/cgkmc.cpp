#!/usr/bin/env node

import fs from "node:fs";
import path from "node:path";

const resultsDir = path.join(process.cwd(), "results");
const logPath = path.join(resultsDir, "petn.log");
const ovitoSurfacePath = path.join(resultsDir, "petn_ovito_surface_area.csv");
const proxyCsvPath = path.join(resultsDir, "petn_surface_energy.csv");
const outputPath = path.join(resultsDir, "petn_surface_energy_ovito.csv");

const totalSites = 62 * 62 * 140 * 2;
const cohesiveEnergyEV = -1.034;
const evPerA2ToMJPerM2 = 16021.76634;
const defaultPetnKappa = 3009322.460;

function argValue(name, fallback) {
  const index = process.argv.indexOf(name);
  return index >= 0 && process.argv[index + 1] ? process.argv[index + 1] : fallback;
}

const kappa = Number(argValue("--kappa", process.env.PETN_KAPPA ?? defaultPetnKappa));
if (!Number.isFinite(kappa) || kappa <= 0) {
  throw new Error(`Invalid kappa value: ${kappa}`);
}

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

function existingTimeKappaByFrame() {
  if (!fs.existsSync(proxyCsvPath)) {
    return new Map();
  }
  return new Map(
    parseCsv(proxyCsvPath).map((row) => [Number(row.frame), Number(row.time_kappa)])
      .filter(([frame, timeKappa]) => Number.isFinite(frame) && Number.isFinite(timeKappa))
  );
}

const logEntries = parseLog();
const surfaceRows = parseCsv(ovitoSurfacePath);
const proxyTimeKappa = existingTimeKappaByFrame();

if (surfaceRows.length !== logEntries.length) {
  throw new Error(`Frame/log mismatch: OVITO CSV has ${surfaceRows.length} rows, log has ${logEntries.length}`);
}

const outputRows = [[
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
  const step = row.timestep && row.timestep !== "undefined" ? row.timestep : frame;
  const timeKappa = proxyTimeKappa.get(frame) ?? log.time * kappa;

  if (!log || !Number.isFinite(surfaceArea) || surfaceArea <= 0 || !Number.isFinite(nAtoms)) {
    throw new Error(`Invalid OVITO surface row for frame ${frame}: ${JSON.stringify(row)}`);
  }

  const gammaEVPerA2 = (log.energy - cohesiveEnergyEV * nAtoms) / surfaceArea;
  outputRows.push([
    frame,
    step,
    log.time,
    timeKappa,
    log.energy,
    nAtoms,
    surfaceArea,
    gammaEVPerA2,
    gammaEVPerA2 * evPerA2ToMJPerM2,
  ].join(","));
}

fs.writeFileSync(outputPath, `${outputRows.join("\n")}\n`);
console.log(`Wrote ${outputPath}`);
