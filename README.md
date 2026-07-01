# CGkMC.cpp

`cgkmc.cpp` is a single-threaded C++20 implementation of the Crystal Growth Kinetic
Monte Carlo (`cgkmc`) model. It stores interactions as compact neighbor
lists and updates local interface state incrementally after each KMC event.

The main architectural difference from the 'cgkmc' Python implementation is that each
KMC step avoids a full sparse-matrix interface rebuild over the lattice. The C++
core maintains:

- occupied sites
- occupied-neighbor counts
- local energy fields
- solid and solvent interface sets
- Fenwick-tree evaporation weights

## Requirements

- CMake 3.20 or newer
- C++20 compiler, such as AppleClang, Clang, or GCC
- Python with the `ovito` package
- Node.js
- Gnuplot

## Build

From the repository root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

This builds:

- `build/libcgkmc_core.a`
- `build/cgkmc_run`

The positional arguments are:

```text
cgkmc_run <input_json> <dump_file> <dump_every> [log_file]
```

The optional log file writes lines with `TIME=... ENERGY=... OCCUPANCY=...`,
matching the downstream figure scripts.

## Run the PETN Simulation

```bash
mkdir -p results

build/cgkmc_run \
  examples/petn_paper.json \
  results/petn.dump \
  1000 \
  results/petn.log
```

The PETN input uses:

- lattice: `62 x 62 x 140` unit cells
- basis: `(0, 0, 0)` and `(1/2, 1/2, 1/2)`
- lattice parameters: `a = b = 9.087 A`, `c = 6.738 A`
- cutoffs: `7.0`, `7.5`, `9.5 A`
- interaction energies: `-0.294`, `-0.184`, `-0.002 eV`
- initial radius: `75 A`
- KMC steps: `1,000,000`

## Run the L-Tryptophan Example

The TRP input is a density-preserving orthorhombic surrogate for the reported
triclinic L-tryptophan crystal. See
`examples/trp_l_tryptophan_notes.md` before interpreting morphology results.

```bash
mkdir -p results

build/cgkmc_run \
  examples/trp_l_tryptophan.json \
  results/trp_l_tryptophan.dump \
  1000 \
  results/trp_l_tryptophan.log
```

After exporting OVITO surface areas for TRP, generate TRP figures with:

```bash
python tools/export_ovito_surface_area.py \
  --dump results/trp_l_tryptophan.dump \
  --output results/trp_l_tryptophan_ovito_surface_area.csv \
  --radius 7.0 \
  --smoothing 100

node tools/build_trp_figures.mjs
gnuplot tools/plot_trp_figure2_ovito.gp
gnuplot tools/plot_trp_figure3.gp
```

Outputs:

- `results/trp_l_tryptophan_surface_energy_ovito.csv`
- `results/figure2_trp_l_tryptophan_surface_energy.png`
- `results/figure2_trp_l_tryptophan_surface_energy_time_s.png`
- `results/trp_l_tryptophan_final_points.tsv`
- `results/trp_l_tryptophan_final_surface_points.tsv`
- `results/figure3_trp_l_tryptophan_final_morphology.png`

## Recreate Figure 2 with OVITO Surface Areas

Run this after generating `results/petn.dump` and `results/petn.log`:

```bash
python -m pip install ovito

python tools/export_ovito_surface_area.py \
  --dump results/petn.dump \
  --output results/petn_ovito_surface_area.csv \
  --radius 7.0 \
  --smoothing 100

node tools/build_petn_figure2_from_ovito.mjs

gnuplot tools/plot_petn_figure2_ovito.gp
```

Outputs:

- `results/petn_ovito_surface_area.csv`
- `results/petn_surface_energy_ovito.csv`
- `results/figure2_petn_surface_energy.png`
- `results/figure2_petn_surface_energy_full.png`

## Recreate Figure 3

Run this after generating `results/petn.dump` and `results/petn.log`:

```bash
node tools/recreate_petn_figures.mjs
gnuplot tools/plot_petn_figure3.gp
```

Outputs:

- `results/petn_surface_energy.csv`
- `results/petn_final_points.tsv`
- `results/petn_final_surface_points.tsv`
- `results/figure3_petn_final_morphology.svg`
- `results/figure3_petn_final_morphology.png`

The Node helper derives the final-frame point cloud and surface-site point cloud
from the C++ dump. The Gnuplot script renders the PNG view.

## Citation and Attribution

This project is a C++ implementation of the model and PETN reproduction workflow
from the original `cgkmc` work. If you use this repository, cite:

```bibtex
@article{jeffries2026petnCgkmc,
  title = {Kinetic Monte Carlo Prediction of the Morphology of Pentaerythritol Tetranitrate},
  author = {Jeffries, Jacob and Singh, Himanshu and Perriot, Romain and Negre, Christian and Redondo, Antonio and Martinez, Enrique},
  journal = {Crystal Growth \& Design},
  year = {2026},
  doi = {10.1021/acs.cgd.5c01428}
}
```

Please also acknowledge the original Python package:

```text
Crystal Growth Kinetic Monte Carlo (cgkmc)
https://github.com/jwjeffr/cgkmc
```

The PETN Figure 2 and Figure 3 reproduction commands in this repository are
intended to validate against the same paper-scale PETN setup.

## License

MIT. See `LICENSE`.
