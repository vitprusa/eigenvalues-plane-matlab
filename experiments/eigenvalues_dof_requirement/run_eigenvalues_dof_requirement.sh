#!/usr/bin/env bash
# How many degrees of freedom the leading eigenvalues cost, with DST, FD and FEM.
#
# Sweeps the resolutions of each method, measures the largest relative error over
# the first n eigenvalues, and reports the dof count needed to hold that error
# under a tolerance. The run table, the answer table and the figure go into
# results/eigenvalues_dof_requirement/.
#
#   ./run_eigenvalues_dof_requirement.sh                 # first 1000, every domain
#   ./run_eigenvalues_dof_requirement.sh 1000 L_shaped   # just that one
#
# A third argument sets the relative accuracies asked for, as a space-separated
# list; the default is 0.1 0.01 0.001.
#
#   ./run_eigenvalues_dof_requirement.sh 1000 L_shaped "0.08 0.04 0.02 0.01 0.001"
#
# Existing CSVs are reused run by run (a figure-only pass costs nothing), so a
# new set of tolerances costs nothing either; delete them to recompute.
# Requires the PDE Toolbox.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
count="${1:-1000}"
name="${2:-}"
tolerances="${3:-}"

if ! command -v matlab >/dev/null 2>&1; then
	echo "Error: matlab not found. Install MATLAB and ensure it is on your PATH." >&2
	exit 1
fi

matlab -batch "addpath('${script_dir}'); compute_eigenvalues_dof_requirement(${count}, '${name}', [${tolerances}])"
