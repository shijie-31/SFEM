function K = cal_K(elem, node)
%% Pre‑processing ---------------------------------------------------------
nDof = size(node, 1);
nEl  = numel(elem);
% Rough estimate of non‑zeros:   ~20 entries per element (tweak if needed)
K = spalloc(nDof, nDof, nEl * 20);

%% Element loop -----------------------------------------------------------
for iEl = 1 : nEl
    % ---- Element data ---------------------------------------------------
    elNodes  = elem{iEl};             % node indices of the current element
    nENodes  = numel(elNodes);        % number of nodes in the element
    elCoords = node(elNodes, :);      % nodal coordinates (nENodes × 2)
    % ---- Fan triangulation of the polygonal element --------------------
    [V, T, nTri] = fan_triangulation(elCoords, 1);
    [midPts, areas, edgeLens, edgeNorms] = tri_info(V, T);
    Ke = sparse(nENodes, nENodes);
    for iTri = 1 : nTri
        C     = squeeze(midPts(iTri, :, :));    % 3 mid‑points of current triangle
        nVec  = squeeze(edgeNorms(iTri, :, :)); % 3 edge normals (unit vectors)
        lEdge = edgeLens(iTri, :);              % 3 edge lengths
        area  = areas(iTri);                    % triangle area
        % Compute RPIM shape functions at the three mid‑points
        N = zeros(3, nENodes);                  % pre‑allocate shape‑function matrix
        for g = 1 : 3
            N(g, :) = RPIM(elCoords, C(g, :));
        end
        % Gradient matrix B (3 × nENodes)
        B = (lEdge(1) * nVec(1, :)' * N(1, :) + ...
             lEdge(2) * nVec(2, :)' * N(2, :) + ...
             lEdge(3) * nVec(3, :)' * N(3, :)) / area;
        Ke = Ke + B' * B * area;
    end
    % ---- Assembly into global stiffness matrix -------------------------
    dofIdx = get_eledof(elNodes, nENodes, 1);   % system DOF indices for the element
    K(dofIdx, dofIdx) = K(dofIdx, dofIdx) + Ke; % add elemental matrix to global
end
end
