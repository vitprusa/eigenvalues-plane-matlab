function [evals, info] = fem_laplace_spectrum(entry, method)
%FEM_LAPLACE_SPECTRUM Dirichlet-Laplacian spectrum on a domain via the PDE Toolbox.
%
%   [evals, info] = fem_laplace_spectrum(entry, method)
%
%   Builds the geometry described by the catalog entry (see DOMAIN_CATALOG_FEM)
%   with decsg, meshes it with quadratic triangles, imposes homogeneous
%   Dirichlet conditions on every edge, and solves -Laplacian u = lambda u by
%   one of two methods.
%
%   Inputs:
%     entry  - struct with fields name, gd, ns, sf, Hmax_eig, Hmax_solvepdeeig
%              (see DOMAIN_CATALOG_FEM).
%     method - "eig"         : assemble the stiffness/mass matrices with
%                              assembleFEMatrices and solve the dense
%                              generalized problem eig(K, M); mesh size
%                              entry.Hmax_eig (coarser, dense solve).
%              "solvepdeeig" : high-level solvepdeeig over the range [0, 200];
%                              mesh size entry.Hmax_solvepdeeig (finer).
%
%   Outputs:
%     evals - eigenvalues of -Laplacian, ascending. The low modes are the
%             accurate ones.
%     info  - struct with fields method, Hmax, dofs, n_nodes, time. For both
%             methods dofs is the number of free (unconstrained) degrees of
%             freedom -- size(K,1) for "eig", size(B,2) of the nullspace basis
%             for "solvepdeeig" -- while n_nodes is the raw mesh-node count.
%             time excludes the extra "solvepdeeig" nullspace assembly.
%
%   Requires the PDE Toolbox.
%
%   See also DOMAIN_CATALOG_FEM, COMPUTE_SPECTRUM_FEM_EIG,
%   COMPUTE_SPECTRUM_FEM_SOLVEPDEEIG.

    method = string(method);
    t0 = tic;

    % Geometry and model set-up (shared by both methods).
    dl = decsg(entry.gd, entry.sf, entry.ns);
    model = createpde();
    geometryFromEdges(model, dl);
    applyBoundaryCondition(model, 'dirichlet', 'Edge', 1:model.Geometry.NumEdges, 'u', 0);
    % -div(c*grad u) + a*u = lambda*d*u, with c = d = 1 and a = 0.
    specifyCoefficients(model, 'm', 0, 'd', 1, 'c', 1, 'a', 0, 'f', 0);

    switch method
        case "eig"
            Hmax = entry.Hmax_eig;
            generateMesh(model, 'Hmax', Hmax, 'GeometricOrder', 'quadratic');
            % Assemble K, M and remove the constrained DOFs with the nullspace
            % basis B, then solve the dense generalized eigenproblem.
            FEM_raw = assembleFEMatrices(model, 'KM');
            FEM_ns  = assembleFEMatrices(model, 'nullspace');
            B = FEM_ns.B;
            K = B' * FEM_raw.K * B;
            Mmat = B' * FEM_raw.M * B;
            evals = sort(real(eig(full(K), full(Mmat))), 'ascend');
            dofs = size(K, 1);
            elapsed = toc(t0);

        case "solvepdeeig"
            Hmax = entry.Hmax_solvepdeeig;
            generateMesh(model, 'Hmax', Hmax, 'GeometricOrder', 'quadratic');
            result = solvepdeeig(model, [0, 200]);
            evals = sort(real(result.Eigenvalues), 'ascend');
            elapsed = toc(t0);   % stop timing before the extra DOF-count assembly
            % True (free) DOF count via the nullspace basis, matching the "eig"
            % path's size(K,1). This assembly is deliberately excluded from the
            % reported time, which covers only the solvepdeeig solve.
            FEM_ns = assembleFEMatrices(model, 'nullspace');
            dofs = size(FEM_ns.B, 2);

        otherwise
            error('fem_laplace_spectrum:badMethod', ...
                'method must be "eig" or "solvepdeeig", got "%s".', method);
    end

    info = struct('method', method, 'Hmax', Hmax, 'dofs', dofs, ...
                  'n_nodes', size(model.Mesh.Nodes, 2), 'time', elapsed);
end
