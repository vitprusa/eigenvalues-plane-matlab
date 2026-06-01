#!/usr/bin/env bash
# Run the Wolfram Language reference eigenvalue computation, writing the
# *-eigenvalues.csv files into results/wolfram/.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v wolframscript >/dev/null 2>&1; then
	echo "Error: wolframscript not found. Install Wolfram Engine or Mathematica and ensure wolframscript is on your PATH." >&2
	exit 1
fi

wolframscript -file "${script_dir}/compute_wolfram_spectra.wls"
