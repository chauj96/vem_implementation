function V3 = perturbVertices(V3, n, L, amp, seed, freezePlanes)
% randomly displace interior vertices of a structured box mesh
%
%   amp          : displacement amplitude as a fraction of the local cell size
%   freezePlanes : nP x 2 list of [dim, value] planes to keep fixed, on top of
%                  the domain boundary. defaults to the material interface
%                  planes x = 1/2, y = 1/2, z = 1/2 of section 5.1.1, so a
%                  perturbed grid still resolves the heterogeneity exactly
%
% vertices on the boundary or on a frozen plane are not moved, and directions
% with a single cell layer stay flat so quasi-2D slabs remain slabs.
%
% NOTE: this moves vertices only. quadrilateral faces of a hexahedral mesh
% become non-planar under perturbation, so cartesianMeshBox splits them into
% triangles; see its 'planar' option.

    if nargin < 6 || isempty(freezePlanes)
        freezePlanes = [1 0.5; 2 0.5; 3 0.5];
    end

    if amp <= 0
        return;
    end

    n = n(:).';
    h = L(:).' ./ max(n, 1);

    rng(seed, 'twister');

    tol = 1e-12;
    movable = true(size(V3,1), 3);
    active = (n > 1);

    for d = 1:3
        onBoundary = abs(V3(:,d)) < tol | abs(V3(:,d) - L(d)) < tol;
        movable(onBoundary, d) = false;
    end

    % keep the material interface exactly on mesh faces
    for p = 1:size(freezePlanes,1)
        d = freezePlanes(p,1);
        v = freezePlanes(p,2);
        if v > tol && v < L(d) - tol
            movable(abs(V3(:,d) - v) < tol, d) = false;
        end
    end

    r = (rand(size(V3)) - 0.5) * 2;      % in [-1,1]

    for d = 1:3
        if active(d)
            m = movable(:,d);
            V3(m,d) = V3(m,d) + amp * h(d) * r(m,d);
        end
    end

end
