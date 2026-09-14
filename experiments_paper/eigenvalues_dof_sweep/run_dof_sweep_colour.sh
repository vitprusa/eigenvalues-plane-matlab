#!/usr/bin/env bash
# Regenerate the colour paper DOF-sweep figures (no titles; method by colour,
# DOF level by line style, thicker lines) from the existing CSVs in
# results/eigenvalues_dof_sweep/, writing EPS figures into
# results_paper/eigenvalues_dof_sweep/. No recompute: the CSVs must already be
# present (run the compute_* scripts in experiments/eigenvalues_dof_sweep first
# if they are not).
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v matlab >/dev/null 2>&1; then
	echo "Error: matlab not found. Install MATLAB and ensure it is on your PATH." >&2
	exit 1
fi

matlab -batch "addpath('${script_dir}'); \
	plot_rectangle_dof_sweep_colour(); \
	plot_isosceles_triangle_dof_sweep_colour(); \
	plot_L_shaped_dof_sweep_colour()"
