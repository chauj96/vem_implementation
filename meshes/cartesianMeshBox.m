function [cell_struct, face_struct, V3, cells3D] = cartesianMeshBox(n, L, perturb, seed, planar)
% structured hexahedral mesh of the box [0,L(1)] x [0,L(2)] x [0,L(3])
%
%   n       : [nx ny nz] cells per direction
%   L       : [Lx Ly Lz]   (default [1 1 1])
%   perturb : relative amplitude of random interior-vertex perturbation,
%             0 for a regular grid (default 0). directions with n == 1 are
%             left unperturbed so quasi-2D slabs stay flat
%   seed    : rng seed, so all methods run on the same perturbed grid
%   planar  : split quadrilateral faces into triangles so every facet is
%             planar. defaults to true when perturb > 0, since perturbing a
%             hexahedral mesh makes its quad faces non-planar. a perturbed hex
%             then becomes a 12-face polyhedron

    if nargin < 2 || isempty(L), L = [1 1 1]; end
    if nargin < 3 || isempty(perturb), perturb = 0; end
    if nargin < 4 || isempty(seed), seed = 0; end
    if nargin < 5 || isempty(planar), planar = (perturb > 0); end

    n = n(:).'; if isscalar(n), n = [n n n]; end
    nx = n(1); ny = n(2); nz = n(3);

    [X, Y, Z] = ndgrid(linspace(0,L(1),nx+1), linspace(0,L(2),ny+1), linspace(0,L(3),nz+1));
    V3 = [X(:), Y(:), Z(:)];

    vid = @(i,j,k) i + (nx+1)*(j-1) + (nx+1)*(ny+1)*(k-1);

    nfx = (nx+1)*ny*nz;
    nfy = nx*(ny+1)*nz;
    nfz = nx*ny*(nz+1);
    nFaces = nfx + nfy + nfz;

    fx = @(i,j,k) i + (nx+1)*(j-1) + (nx+1)*ny*(k-1);
    fy = @(i,j,k) nfx + i + nx*(j-1) + nx*(ny+1)*(k-1);
    fz = @(i,j,k) nfx + nfy + i + nx*(j-1) + nx*ny*(k-1);

    cid = @(i,j,k) i + nx*(j-1) + nx*ny*(k-1);

    faceVerts = cell(nFaces,1);
    faceCells = zeros(nFaces,2);

    % x-normal faces
    for k = 1:nz, for j = 1:ny, for i = 1:nx+1
        f = fx(i,j,k);
        faceVerts{f} = [vid(i,j,k), vid(i,j+1,k), vid(i,j+1,k+1), vid(i,j,k+1)];
        faceCells(f,:) = [pick(i>1, cid(max(i-1,1),j,k)), pick(i<=nx, cid(min(i,nx),j,k))];
    end, end, end

    % y-normal faces
    for k = 1:nz, for j = 1:ny+1, for i = 1:nx
        f = fy(i,j,k);
        faceVerts{f} = [vid(i,j,k), vid(i,j,k+1), vid(i+1,j,k+1), vid(i+1,j,k)];
        faceCells(f,:) = [pick(j>1, cid(i,max(j-1,1),k)), pick(j<=ny, cid(i,min(j,ny),k))];
    end, end, end

    % z-normal faces
    for k = 1:nz+1, for j = 1:ny, for i = 1:nx
        f = fz(i,j,k);
        faceVerts{f} = [vid(i,j,k), vid(i+1,j,k), vid(i+1,j+1,k), vid(i,j+1,k)];
        faceCells(f,:) = [pick(k>1, cid(i,j,max(k-1,1))), pick(k<=nz, cid(i,j,min(k,nz)))];
    end, end, end

    cellFaces = cell(nx*ny*nz,1);
    for k = 1:nz, for j = 1:ny, for i = 1:nx
        cellFaces{cid(i,j,k)} = [fx(i,j,k), fx(i+1,j,k), ...
                                 fy(i,j,k), fy(i,j+1,k), ...
                                 fz(i,j,k), fz(i,j,k+1)];
    end, end, end

    V3 = perturbVertices(V3, n, L, perturb, seed);

    if planar
        [faceVerts, faceCells, cellFaces] = triangulateFaces(faceVerts, faceCells, cellFaces);
    end

    [cell_struct, face_struct, cells3D] = buildMeshStructs(V3, faceVerts, faceCells, cellFaces);

end

function v = pick(cond, val)
    if cond, v = val; else, v = 0; end
end
