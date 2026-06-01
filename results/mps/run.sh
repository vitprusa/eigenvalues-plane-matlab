#!/usr/bin/env bash
#
# run_Ldrum_MPS.sh - run Ldrum_modified.m and export the eigenvalues to
# L_shaped_eigenvalues_MPS.csv, using the same CSV format as the files
# produced by the Wolfram Language script (header "n","lambda_n").

set -euo pipefail

# Work from the directory containing this script (where the .m files live).
cd "$(dirname "$0")"

OUT="L_shaped_eigenvalues_MPS.csv"

matlab -nodisplay -nosplash -batch "Ldrum_modified; fid = fopen('${OUT}', 'w'); fprintf(fid, '\"n\",\"lambda_n\"\n'); for n = 1:numel(evals), fprintf(fid, '%d,%.15g\n', n, evals(n)); end; fclose(fid); disp(['Wrote ', num2str(numel(evals)), ' eigenvalues to ${OUT}']);"
