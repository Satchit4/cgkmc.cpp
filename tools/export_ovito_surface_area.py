#!/usr/bin/env python3
"""Export OVITO alpha-shape surface areas for the PETN trajectory."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


SURFACE_AREA_ATTRIBUTES = (
    "ConstructSurfaceMesh.surface_area",
    "ConstructSurfaceMesh.surface_area_total",
)
FILLED_VOLUME_ATTRIBUTES = (
    "ConstructSurfaceMesh.filled_volume",
    "ConstructSurfaceMesh.solid_volume",
)


def first_attribute(attributes, names):
    for name in names:
        if name in attributes:
            return attributes[name]
    return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export per-frame OVITO Construct Surface Mesh areas for cgkmc PETN dumps."
    )
    parser.add_argument("--dump", type=Path, default=Path("results/petn.dump"))
    parser.add_argument("--output", type=Path, default=Path("results/petn_ovito_surface_area.csv"))
    parser.add_argument("--radius", type=float, default=7.0)
    parser.add_argument("--smoothing", type=int, default=100)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    try:
        from ovito.io import import_file
        from ovito.modifiers import ConstructSurfaceModifier
    except ImportError as exc:
        raise SystemExit(
            "OVITO's Python module is required. Install it in the active environment with:\n"
            "  python -m pip install ovito\n"
            "or run this file with OVITO's bundled ovitos interpreter."
        ) from exc

    pipeline = import_file(str(args.dump))
    pipeline.modifiers.append(
        ConstructSurfaceModifier(
            method=ConstructSurfaceModifier.Method.AlphaShape,
            radius=args.radius,
            smoothing_level=args.smoothing,
            identify_regions=True,
        )
    )

    num_frames = getattr(pipeline, "num_frames", None)
    if num_frames is None:
        num_frames = pipeline.source.num_frames

    with args.output.open("w", newline="") as output_file:
        writer = csv.writer(output_file)
        writer.writerow([
            "frame",
            "timestep",
            "n_atoms",
            "surface_area_A2",
            "filled_volume_A3",
        ])

        for frame in range(num_frames):
            data = pipeline.compute(frame)
            surface_area = first_attribute(data.attributes, SURFACE_AREA_ATTRIBUTES)
            if surface_area is None:
                available = ", ".join(sorted(data.attributes.keys()))
                raise RuntimeError(
                    "OVITO did not report a recognized ConstructSurfaceMesh surface-area attribute. "
                    f"Available attributes: {available}"
                )

            particle_count = getattr(data.particles, "count", None)
            if particle_count is None:
                particle_count = len(data.particles.positions)

            writer.writerow([
                frame,
                data.attributes.get("Timestep", frame),
                particle_count,
                surface_area,
                first_attribute(data.attributes, FILLED_VOLUME_ATTRIBUTES) or "",
            ])

            if frame == 0 or (frame + 1) % 50 == 0 or frame + 1 == num_frames:
                print(f"Processed OVITO surface frame {frame + 1} / {num_frames}")

    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
