function [G, fractureInfo] = inclinedFracture(Lx, Ly, Lz, h, Nz, psi, b, Nfrac)
% fracture-conforming graded triangular prism mesh
%
% domain: [-Lx/2,Lx/2] x [-Ly/2,Ly/2] x [0,Lz]
%
% fracture: straight line centered at origin half-length b, inclination psi (degrees)

    %% 1. Geometry
    
    assert(all([Lx, Ly, Lz, h, b] > 0), 'All geometric parameters must be positive.');
    
    assert(Nz >= 1 && Nz == round(Nz), 'Nz must be a positive integer.');
    
    theta = deg2rad(psi);
    
    tangent2D = [cos(theta), sin(theta)];
    normal2D = [-sin(theta), cos(theta)];
    
    tangent = [tangent2D, 0]';
    normal = [normal2D, 0]';
    
    p1 = -b*tangent2D;
    p2 = b*tangent2D;
    
    assert(all(abs([p1; p2]) < [Lx/2, Ly/2], 'all'), 'Fracture endpoints must lie inside domain.');
    
    hFrac = 2 * b/Nfrac;
    
    % local and far-field target spacing
    % hNear = 0.85*hFrac;
    % hFar = max(h, 1.5*hNear);
    
    tol = 1e-8 * max([Lx, Ly, Lz, 1]);
    
    
    %% 2. Fracture vertices
    
    s = linspace(-b, b, Nfrac+1)';
    Pfrac = s * tangent2D;
    
    %% 3. Generate more regular graded background points
    
    rng(11);
    
    hNear = 0.85 * hFrac;
    hFar = max(h, 1.5*hNear);
    
    % generate candidates on a triangular lattice
    step = hNear;
    dy = sqrt(3)/2*step;
    
    ys = (-Ly/2+step):dy:(Ly/2-step);
    
    Pint = zeros(0,2);
    
    for j = 1:numel(ys)
    
        y = ys(j);
    
        shift = mod(j,2)*step/2;
    
        xs = (-Lx/2+step+shift):step:(Lx/2-step);
    
        X = [xs(:), y*ones(numel(xs),1)];
    
        d = distanceToSegment(X,p1,p2);
    
        % hTarget = hNear + (hFar-hNear)*(1-exp(-d/(2.5*b)));
        refineWidth = 0.5*b;
    
        hTarget = hNear + (hFar-hNear)*(1-exp(-(d/refineWidth).^2));
    
        % deterministic spatial filtering
        for i = 1:size(X,1)
    
            candidate = X(i,:);
    
            if isempty(Pint)
                Pint = candidate;
                continue;
            end
    
            distExisting = sqrt(sum((Pint-candidate).^2,2));
    
            if all(distExisting > 0.85*hTarget(i))
                Pint(end+1,:) = candidate; 
            end
        end
    end
    
    %% 4. Remove points too close to fracture
    
    d = distanceToSegment(Pint,p1,p2);
    Pint = Pint(d > 0.55*hNear,:);
    
    
    %% 5. Generate rectangular boundary points
    
    nx = max(2,ceil(Lx/hFar));
    ny = max(2,ceil(Ly/hFar));
    
    xv = linspace(-Lx/2,Lx/2,nx+1)';
    yv = linspace(-Ly/2,Ly/2,ny+1)';
    
    bottom = [xv, -Ly/2*ones(numel(xv),1)];
    
    right = [Lx/2*ones(ny-1,1), yv(2:end-1)];
    
    top = [flipud(xv), Ly/2*ones(numel(xv),1)];
    
    left = [-Lx/2*ones(ny-1,1), flipud(yv(2:end-1))];
    
    Pboundary = [bottom; right; top; left];
    
    nb = size(Pboundary,1);
    
    boundaryEdges = [(1:nb)', [(2:nb)'; 1]];
    
    
    %% 6. Combine points and constrained edges
    
    P = [Pboundary; Pfrac; Pint];
    
    fracIds = (nb+1:nb+Nfrac+1)';
    fractureEdges = [fracIds(1:end-1), fracIds(2:end)];
    constraints = [boundaryEdges; fractureEdges];
    
    
    %% 7. Constrained Delaunay triangulation
    
    DT = delaunayTriangulation(P,constraints);
    
    P = DT.Points;
    T = DT.ConnectivityList;
    
    % verify that all fracture segments remain mesh edges
    triEdges = [ ...
        T(:,[1,2]); ...
        T(:,[2,3]); ...
        T(:,[3,1])];
    
    triEdges = sort(triEdges,2);
    
    [uniqueEdges,~,ic] = unique(triEdges,'rows');
    
    edgeCounts = accumarray(ic,1);
    
    requiredEdges = sort(fractureEdges,2);
    
    [found,locations] = ismember(requiredEdges,uniqueEdges,'rows');
    
    assert(all(found), 'Some constrained fracture edges are missing.');
    assert(all(edgeCounts(locations)==2), 'Some fracture edges are not internal edges.');
    
    
    %% 8. Convert triangulation to MRST grid
    
    G2D = triangleGrid(P,T);
    G2D = computeGeometry(G2D);
    
    assert(G2D.cells.num == size(T,1), 'MRST triangle conversion failed.');
    
    
    %% 9. Identify 2D fracture faces
    
    faces2D = zeros(0,1);
    
    for f = 1:G2D.faces.num
    
        % only internal edges
        if any(G2D.faces.neighbors(f,:)==0)
            continue;
        end
    
        ids = G2D.faces.nodes(G2D.faces.nodePos(f): G2D.faces.nodePos(f+1)-1);
    
        xy = G2D.nodes.coords(ids,1:2);
    
        if size(xy,1)~=2
            continue;
        end
    
        % projection onto fracture tangent and normal
        along = xy*tangent2D';
        perp  = xy*normal2D';
    
        onLine = all(abs(perp)<tol);
    
        onSegment = all(along>=-b-tol) && all(along<=b+tol);
    
        if onLine && onSegment
    
            faces2D(end+1,1) = f;
    
        end
    end
    
    assert(numel(faces2D)==Nfrac, 'Expected %d fracture edges, found %d.', Nfrac,numel(faces2D));
    
    fractureLength = sum(G2D.faces.areas(faces2D));
    
    assert(abs(fractureLength-2*b)<1e-7*max(1,2*b), 'Fracture length mismatch.');
    
    
    %% 10. Extrude triangular mesh into 3D
    
    G = makeLayeredGrid(G2D,Nz);
    
    z = G.nodes.coords(:,3);
    
    zmin = min(z);
    zmax = max(z);
    
    assert(zmax>zmin, 'Invalid extrusion thickness.');
    
    G.nodes.coords(:,3) = Lz*(z-zmin)/(zmax-zmin);
    
    G = computeGeometry(G);
    
    assert(G.cells.num==G2D.cells.num*Nz, 'Unexpected number of 3D cells.');
    
    
    %% 11. Identify 3D fracture faces
    
    fractureFaces = zeros(0,1);
    
    for f = 1:G.faces.num
    
        % fracture faces must be internal
        if any(G.faces.neighbors(f,:)==0)
            continue;
        end
    
        ids = G.faces.nodes(G.faces.nodePos(f): G.faces.nodePos(f+1)-1);
    
        xyz = G.nodes.coords(ids,:);
    
        % vertical prism faces have four vertices
        if size(xyz,1)~=4
            continue;
        end
    
        % skip horizontal faces
        if max(xyz(:,3))-min(xyz(:,3))<tol
            continue;
        end
    
        xy = xyz(:,1:2);
    
        along = xy*tangent2D';
        perp = xy*normal2D';
    
        onLine = all(abs(perp)<tol);
    
        onSegment = all(along>=-b-tol) && all(along<=b+tol);
    
        if onLine && onSegment
    
            fractureFaces(end+1,1) = f; 
    
        end
    end
    
    expectedFaces = Nfrac*Nz;
    
    assert(numel(fractureFaces)==expectedFaces, 'Expected %d 3D fracture faces, found %d.', expectedFaces,numel(fractureFaces));
    
    
    %% 12. Identify plus/minus neighboring cells
    
    fractureCells = G.faces.neighbors(fractureFaces,:);
    
    plusCells = zeros(numel(fractureFaces),1);
    minusCells = zeros(numel(fractureFaces),1);
    
    for k = 1:numel(fractureFaces)
    
        c1 = fractureCells(k,1);
        c2 = fractureCells(k,2);
    
        x1 = G.cells.centroids(c1,1:2);
        x2 = G.cells.centroids(c2,1:2);
    
        d1 = dot(x1,normal2D);
        d2 = dot(x2,normal2D);
    
        assert(abs(d1-d2)>tol, 'Ambiguous fracture orientation.');
    
        if d1>d2
    
            plusCells(k) = c1;
            minusCells(k) = c2;
    
        else
    
            plusCells(k) = c2;
            minusCells(k) = c1;
    
        end
    end
    
    
    %% 13. Store fracture information
    
    fractureInfo.center = [0;0;0];
    
    fractureInfo.tangent = tangent;
    fractureInfo.normal = normal;
    
    fractureInfo.psi = psi;
    fractureInfo.b = b;
    
    fractureInfo.endpoints = [p1;p2];
    
    fractureInfo.faces2D = faces2D;
    fractureInfo.faces = fractureFaces;
    
    fractureInfo.nFaces = numel(fractureFaces);
    
    fractureInfo.cells = fractureCells;
    fractureInfo.plusCells = plusCells;
    fractureInfo.minusCells = minusCells;
    
    fractureInfo.area = sum(G.faces.areas(fractureFaces));
    
    fractureInfo.nSegments = Nfrac;
    
    assert(abs(fractureInfo.area-2*b*Lz) < 1e-7*max(1,2*b*Lz), 'Fracture area mismatch.');
    
    
    %% 14. Summary
    
    fprintf('\n[# TRIANGULAR PRISM FRACTURE MESH]\n');
    
    fprintf('  2D triangles: %d\n',G2D.cells.num);
    fprintf('  3D prisms: %d\n',G.cells.num);
    
    fprintf('  2D fracture edges: %d\n',numel(faces2D));
    fprintf('  3D fracture faces: %d\n',numel(fractureFaces));
    
    fprintf('  Fracture length: %.6f\n',fractureLength);
    fprintf('  Fracture area: %.6f\n',fractureInfo.area);
    
    
    %% 15. Visualize 2D mesh
    
    figure;
    
    plotGrid(G2D,'FaceColor','none');
    
    hold on;
    
    plot([p1(1),p2(1)], [p1(2),p2(2)], 'r-','LineWidth',2);
    
    axis equal;
    
    xlabel('x');
    ylabel('y');
    
    title(sprintf('%d triangles, %d fracture edges', G2D.cells.num,Nfrac));

end

%% HELPER FUNCTION
function d = distanceToSegment(P,a,b)
% distance from points P to finite line segment [a,b].

    a = a(:)';
    b = b(:)';
    
    v = b-a;
    
    tau = ((P-a)*v')/dot(v,v);
    
    tau = max(0,min(1,tau));
    
    closest = a+tau*v;
    
    d = sqrt(sum((P-closest).^2,2));

end

