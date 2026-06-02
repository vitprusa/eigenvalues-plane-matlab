function [evals, dofs] = read_eigs_csv(csv_path)
%READ_EIGS_CSV Read an eigenvalue CSV: the lambda_n column and the DOF count.
%
%   [evals, dofs] = read_eigs_csv(csv_path) parses the "#"-prefixed header for
%   "dofs = <n>" and reads the n, lambda_n table. Returns evals (column vector)
%   and dofs (scalar, NaN if the header has no dofs field).

    dofs = NaN;
    fid = fopen(csv_path, 'r');
    if fid == -1
        error('read_eigs_csv:cannotOpen', 'Could not open %s.', csv_path);
    end
    line = fgetl(fid);
    while ischar(line) && ~isempty(line) && line(1) == '#'
        tok = regexp(line, 'dofs\s*=\s*(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            dofs = str2double(tok{1});
        end
        line = fgetl(fid);
    end
    fclose(fid);

    T = readtable(csv_path, 'CommentStyle', '#');
    evals = T.lambda_n;
end
