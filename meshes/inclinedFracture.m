function [G, fractureInfo] = inclinedFracture(Lx, Ly, Lz, h, Nz, psi, b, Nfrac, meshType)
% fracture-conforming extruded mesh with tetra or hexa options

% tetra: triangular prisms; hexa: composite pebi polygonal prisms

% domain: [-lx/2,lx/2] x [-ly/2,ly/2] x [0,lz]

    if nargin < 9
        meshType = 'hexa';
    end
    meshType = validatestring(meshType, {'tetra', 'hexa'});

    %% 1. geometry

    assert(all([Lx, Ly, Lz, h, b] > 0), 'All geometric parameters must be positive.');

    assert(Nz >= 1 && Nz == round(Nz), 'Nz must be a positive integer.');

    assert(Nfrac >= 1 && Nfrac == round(Nfrac), 'Nfrac must be a positive integer.');

    theta = deg2rad(psi);

    tangent2D = [cos(theta), sin(theta)];
    normal2D = [-sin(theta), cos(theta)];

    tangent = [tangent2D, 0]';
    normal = [normal2D, 0]';

    p1 = -b*tangent2D;
    p2 = b*tangent2D;

    assert(all(abs([p1; p2]) < [Lx/2, Ly/2], 'all'), 'Fracture endpoints must lie inside domain.');

    hFrac = 2*b/Nfrac;

    tol = 1e-8 * max([Lx, Ly, Lz, 1]);

    %% 2. generate 2d mesh

    switch meshType
        case 'tetra'

            %% 2. fracture vertices

            s = linspace(-b, b, Nfrac+1)';
            Pfrac = s * tangent2D;

            %% 3. generate more regular graded background points

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

            %% 4. remove points too close to fracture

            d = distanceToSegment(Pint,p1,p2);
            Pint = Pint(d > 0.55*hNear,:);

            %% 5. generate rectangular boundary points

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

            %% 6. combine points and constrained edges

            P = [Pboundary; Pfrac; Pint];

            fracIds = (nb+1:nb+Nfrac+1)';
            fractureEdges = [fracIds(1:end-1), fracIds(2:end)];
            constraints = [boundaryEdges; fractureEdges];

            %% 7. constrained delaunay triangulation

            DT = delaunayTriangulation(P,constraints);

            P = DT.Points;
            T = DT.ConnectivityList;

            % verify that all fracture segments remain mesh edges
            triEdges = [ T(:,[1,2]); T(:,[2,3]); T(:,[3,1])];

            triEdges = sort(triEdges,2);

            [uniqueEdges,~,ic] = unique(triEdges,'rows');

            edgeCounts = accumarray(ic,1);

            requiredEdges = sort(fractureEdges,2);

            [found,locations] = ismember(requiredEdges,uniqueEdges,'rows');

            assert(all(found), 'Some constrained fracture edges are missing.');
            assert(all(edgeCounts(locations)==2), 'Some fracture edges are not internal edges.');

            %% 8. convert triangulation to mrst grid

            G2D = triangleGrid(P,T);
            G2D = computeGeometry(G2D);

            assert(G2D.cells.num == size(T,1), 'MRST triangle conversion failed.');

            %% 9. identify 2d fracture faces

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
                perp = xy*normal2D';

                onLine = all(abs(perp)<tol);

                onSegment = all(along>=-b-tol) && all(along<=b+tol);

                if onLine && onSegment

                    faces2D(end+1,1) = f;

                end
            end

            assert(numel(faces2D)==Nfrac, 'Expected %d fracture edges, found %d.', Nfrac,numel(faces2D));

            cellLabel = 'triangles';
            prismLabel = 'triangular prisms';

        case 'hexa'

            %% 2. define fracture geometry

            % composite pebi grid uses positive coordinates:
            % [0,lx] x [0,ly]

            shift = [Lx/2, Ly/2];

            % the fracture is specified by its endpoints.
            % the grid generator determines the actual subdivision.

            Pfrac = [p1; p2] + shift;

            faceConstraints = {Pfrac};

            %% 3. generate composite pebi mesh

            % structured background cells with polygonal
            % transition cells near the fracture.

            rng(11);

            G2D = compositePebiGrid2D( [h, h], [Lx, Ly], 'faceConstraints', faceConstraints, 'FCFactor', hFrac/h);

            % return to centered domain

            G2D.nodes.coords(:,1:2) = G2D.nodes.coords(:,1:2) - shift;

            G2D = computeGeometry(G2D);

            assert(all(G2D.cells.volumes > 0), 'Invalid 2D cell areas.');

            %% 4. identify 2d fracture faces

            assert(isfield(G2D.faces,'tag'), 'PEBI grid does not contain face tags.');

            faces2D = find(G2D.faces.tag ~= 0);

            % keep only internal faces
            faces2D = faces2D( all(G2D.faces.neighbors(faces2D,:) > 0,2));

            assert(~isempty(faces2D), 'No tagged fracture edges found.');

            %% snap fracture vertices onto the exact fracture line

            fracNodeIds = [];

            for k = 1:numel(faces2D)

                f = faces2D(k);

                ids = G2D.faces.nodes( G2D.faces.nodePos(f): G2D.faces.nodePos(f+1)-1);

                fracNodeIds = [fracNodeIds; ids(:)];

            end

            fracNodeIds = unique(fracNodeIds);

            % original fracture vertex coordinates
            xy = G2D.nodes.coords(fracNodeIds,1:2);

            % coordinates along and perpendicular to the fracture
            sNode = xy*tangent2D';
            dNode = xy*normal2D';

            fprintf('\n[FRACTURE GEOMETRY CORRECTION]\n');

            fprintf('  Maximum deviation before: %.10e\n', max(abs(dNode)));

            % ensure that only slightly perturbed vertices are corrected
            snapTol = 1e-6*max([Lx,Ly,1]);

            assert(max(abs(dNode)) < snapTol, 'Fracture vertices are too far from the exact line.');

            % identify fracture tip vertices
            [~, iLeft] = min(sNode);
            [~, iRight] = max(sNode);

            assert(abs(sNode(iLeft)+b) < snapTol && abs(sNode(iRight)-b) < snapTol, 'Fracture endpoints mismatch.');

            % set exact fracture endpoints
            sNode(iLeft) = -b;
            sNode(iRight) = b;

            % orthogonal projection onto the exact fracture line
            xyProjected = sNode*tangent2D;

            % update shared mesh vertices
            G2D.nodes.coords(fracNodeIds,1:2) = xyProjected;

            % recompute all geometric quantities
            G2D = computeGeometry(G2D);

            assert(all(G2D.cells.volumes > 0), 'Invalid cell areas after fracture correction.');

            assert(all(G2D.faces.areas > 0), 'Invalid edge lengths after fracture correction.');

            % verify corrected geometry
            dCorrected = G2D.nodes.coords(fracNodeIds,1:2)*normal2D';

            fprintf('  Maximum deviation after:  %.10e\n', max(abs(dCorrected)));

            fprintf('  Maximum node displacement: %.10e\n', max(vecnorm(xyProjected-xy,2,2)));

            cellLabel = 'polygons';
            prismLabel = 'polyhedral prisms';
    end

    %% 3. validate 2d fracture geometry

    nSegments = numel(faces2D);

    fractureLength = sum(G2D.faces.areas(faces2D));

    %% check geometric accuracy

    maxDeviation = 0;
    sMin = inf;
    sMax = -inf;

    for k = 1:nSegments

        f = faces2D(k);

        ids = G2D.faces.nodes( G2D.faces.nodePos(f): G2D.faces.nodePos(f+1)-1);

        xy = G2D.nodes.coords(ids,1:2);

        along = xy*tangent2D';
        perp = xy*normal2D';

        maxDeviation = max(maxDeviation,max(abs(perp)));

        sMin = min(sMin,min(along));
        sMax = max(sMax,max(along));

    end

    fprintf('\n[2D FRACTURE CHECK]\n');

    fprintf('  Target segments: %d\n',Nfrac);
    fprintf('  Tagged segments: %d\n',nSegments);

    fprintf('  Expected length: %.10f\n',2*b);
    fprintf('  Tagged length: %.10f\n',fractureLength);

    fprintf('  Maximum line deviation: %.10e\n',maxDeviation);

    fprintf('  Expected endpoints: [%.6f, %.6f]\n',-b,b);
    fprintf('  Detected endpoints: [%.6f, %.6f]\n',sMin,sMax);

    %% validate fracture geometry

    geomTol = 1e-6*max([Lx,Ly,1]);

    assert(maxDeviation < geomTol, 'Tagged fracture edges deviate from fracture line.');

    assert(abs(sMin+b)<geomTol && abs(sMax-b)<geomTol, 'Fracture endpoints mismatch.');

    assert(abs(fractureLength-2*b)<geomTol, 'Fracture length mismatch.');

    assert(all(G2D.cells.volumes > 0), 'Invalid 2D cell areas.');

    %% 4. extrude mesh into 3d

    G = makeLayeredGrid(G2D,Nz);

    z = G.nodes.coords(:,3);

    zmin = min(z);
    zmax = max(z);

    assert(zmax > zmin, 'Invalid extrusion thickness.');

    G.nodes.coords(:,3) = Lz*(z-zmin)/(zmax-zmin);

    G = computeGeometry(G);

    assert(G.cells.num == G2D.cells.num*Nz, 'Unexpected number of 3D cells.');

    assert(all(G.cells.volumes > 0), 'Invalid 3D cell volumes.');

    %% 5. identify 3d fracture faces

    % match extruded vertical faces to tagged 2d fracture edges.

    fractureFaces = zeros(0,1);

    matchTol = 1e-8 * max([Lx,Ly,Lz,1]);

    % store endpoint coordinates of tagged 2d fracture edges

    edgeXY = cell(nSegments,1);

    for k = 1:nSegments

        f2 = faces2D(k);

        ids = G2D.faces.nodes( G2D.faces.nodePos(f2): G2D.faces.nodePos(f2+1)-1);

        edgeXY{k} = G2D.nodes.coords(ids,1:2);

    end

    % match 3d vertical faces to 2d fracture edges

    for f = 1:G.faces.num

        % only internal faces

        if any(G.faces.neighbors(f,:) == 0)
            continue;
        end

        ids = G.faces.nodes( G.faces.nodePos(f): G.faces.nodePos(f+1)-1);

        xyz = G.nodes.coords(ids,:);

        % skip horizontal faces

        if max(xyz(:,3))-min(xyz(:,3)) < matchTol
            continue;
        end

        % project 3d face vertices onto the xy-plane

        xy = unique(xyz(:,1:2),'rows');

        if size(xy,1) ~= 2
            continue;
        end

        % compare with tagged 2d fracture edges

        for k = 1:nSegments

            edge = edgeXY{k};

            % endpoint matching independent of ordering

            err1 = max(vecnorm(xy-edge,2,2));

            err2 = max(vecnorm(xy-flipud(edge),2,2));

            if min(err1,err2) < matchTol

                fractureFaces(end+1,1) = f;

                break;

            end

        end

    end

    %% validate extrusion of fracture

    expectedFaces = nSegments*Nz;

    fprintf('\n[3D FRACTURE CHECK]\n');

    fprintf('  Expected fracture faces: %d\n',expectedFaces);

    fprintf('  Detected fracture faces: %d\n', numel(fractureFaces));

    assert(numel(fractureFaces) == expectedFaces, 'Expected %d 3D fracture faces, found %d.', expectedFaces,numel(fractureFaces));

    fractureArea = sum(G.faces.areas(fractureFaces));

    expectedArea = fractureLength*Lz;

    fprintf('  Expected fracture area: %.10f\n',expectedArea);

    fprintf('  Detected fracture area: %.10f\n',fractureArea);

    assert(abs(fractureArea-expectedArea) < 1e-7*max(1,expectedArea), 'Extruded fracture area mismatch.');

    %% 6. identify plus/minus neighboring cells

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

        % neighboring cells must lie on opposite
        % sides of the fracture.

        assert((d1 > tol && d2 < -tol) || (d2 > tol && d1 < -tol), 'Ambiguous fracture orientation.');

        if d1 > d2

            plusCells(k) = c1;
            minusCells(k) = c2;

        else

            plusCells(k) = c2;
            minusCells(k) = c1;

        end

    end

    %% 7. store fracture information

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

    fractureInfo.nSegments = nSegments;

    fractureInfo.targetNSegments = Nfrac;

    assert(abs(fractureInfo.area-2*b*Lz) < 1e-7*max(1,2*b*Lz), 'Fracture area mismatch.');

    fractureInfo.meshType = meshType;

    %% 8. summary

    nSides = diff(G2D.cells.facePos);

    fprintf('\n[# %s FRACTURE MESH]\n',upper(meshType));

    fprintf('  2D cells: %d\n',G2D.cells.num);
    fprintf('  3D cells: %d\n',G.cells.num);

    fprintf('  2D fracture edges: %d\n',nSegments);
    fprintf('  3D fracture faces: %d\n',numel(fractureFaces));

    fprintf('  Fracture length: %.6f\n',fractureLength);
    fprintf('  Fracture area: %.6f\n',fractureInfo.area);

    fprintf('  Minimum polygon sides: %d\n',min(nSides));
    fprintf('  Maximum polygon sides: %d\n',max(nSides));
    fprintf('  Average polygon sides: %.2f\n',mean(nSides));

    %% 9. visualize 2d mesh

    figure('Color','w');

    plotGrid(G2D, 'FaceColor','none', 'EdgeColor',[0.4,0.4,0.4]);

    hold on;

    % target fracture geometry
    plot([p1(1),p2(1)], [p1(2),p2(2)], 'r-', 'LineWidth',2);

    % identified fracture edges
    plotFaces(G2D,faces2D, 'EdgeColor','b', 'LineWidth',2);

    axis equal tight;

    xlabel('x');
    ylabel('y');

    title(sprintf('%d %s, %d fracture edges', G2D.cells.num,cellLabel,nSegments));

    %% 10. visualize 3d mesh

    figure('Color','w');

    plotGrid(G, 'FaceColor','none', 'EdgeColor',[0.5,0.5,0.5]);

    hold on;

    plotFaces(G,fractureFaces, 'FaceColor','r', 'FaceAlpha',0.6, 'EdgeColor','r');

    axis equal tight;

    view(35,25);

    xlabel('x');
    ylabel('y');
    zlabel('z');

    title(sprintf('%d %s',G.cells.num,prismLabel));

end

%% helper function
function d = distanceToSegment(P,a,b)
% distance from points p to finite line segment [a,b].

    a = a(:)';
    b = b(:)';

    v = b-a;

    tau = ((P-a)*v')/dot(v,v);

    tau = max(0,min(1,tau));

    closest = a+tau*v;

    d = sqrt(sum((P-closest).^2,2));

end

