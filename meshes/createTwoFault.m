function G = createTwoFault()
    mrstModule add upr

    %% Set boundary
    tx = 1;
    ty = 1;
    tz = 0.4;

    bdr = [0,  0,  0;
           tx, 0,  0;
           tx, ty, 0;
           0,  ty, 0;
           0,  0,  tz;
           tx, 0,  tz;
           tx, ty, tz;
           0,  ty, tz];

    %% Set gridding parameters
    fGs = tx / 16;
    dt  = tx / 40;

    zTilt  = 0.3;
    yBend  = 0.1;
    sinAmp = 0.09;

    rho = @(p) fGs * (1 + 0*p(:,1));

    %% Triangulate the two curved surfaces
    swap1 = [0,1;1,0] == 1;
    swap2 = [1,0;0,1] == 1;

    rectangle1 = [min(bdr(:,[1,3])); max(bdr(:,[1,3]))];

    fixedPts = [rectangle1; ...
                rectangle1(swap1)'; ...
                rectangle1(swap2)'];

    hd = @(p) rho(p)/fGs;
    fd = @(p) drectangle( ...
        p, rectangle1(1), rectangle1(2), ...
        rectangle1(3), rectangle1(4));

    [Pts,t] = distmesh2d(fd, hd, fGs, rectangle1, fixedPts, false);

    t1.ConnectivityList = t;
    t2.ConnectivityList = t;

    faultHeight1z = @(p) ...
        ty/6*ones(size(p,1),1) + zTilt*p(:,2);

    faultHeight2z = @(p) ...
        5*ty/6*ones(size(p,1),1) - zTilt*p(:,2);

    t1.Points = [Pts(:,1), faultHeight1z(Pts), Pts(:,2)];
    t2.Points = [Pts(:,1), faultHeight2z(Pts), Pts(:,2)];

    faultHeight1x = @(p) ...
        p(:,2) + yBend*(p(:,1)-0.75).^2 ...
        + sinAmp*sin(2*pi/tx*p(:,1));

    faultHeight2x = @(p) ...
        p(:,2) + yBend*(p(:,1)-0.75).^2;

    t1.Points(:,2) = faultHeight1x(t1.Points);
    t2.Points(:,2) = faultHeight2x(t2.Points);

    %% Generate conformal surface sites
    gamma = 0.9;

    R1 = circumradiusPerVertex( ...
        t1.Points, t1.ConnectivityList, gamma);

    R2 = circumradiusPerVertex( ...
        t2.Points, t2.ConnectivityList, gamma);

    rho1 = @(p) R1;
    rho2 = @(p) R2;

    F = surfaceSites3D({t1,t2}, {rho1,rho2});

    %% Remove invalid surface sites
    isValid = all(isreal(F.f.pts),2) ...
           & all(isfinite(F.f.pts),2);

    F.f.pts = F.f.pts(isValid,:);
    F.f.Gs  = F.f.Gs(isValid);
    F.f.pri = F.f.pri(isValid);
    F.f.l   = F.f.l(isValid);

    %% Create background reservoir sites
    xmax = max(bdr(:,1)) - dt/2;
    xmin = min(bdr(:,1)) + dt/2;

    ymax = max(bdr(:,2)) - dt/2;
    ymin = min(bdr(:,2)) + dt/2;

    zmax = max(bdr(:,3)) - dt/2;
    zmin = min(bdr(:,3)) + dt/2;

    xr = xmin:dt:xmax;
    yr = ymin:dt:ymax;
    zr = zmin:dt:zmax;

    [X,Y,Z] = ndgrid(xr,yr,zr);

    rSites = [X(:),Y(:),Z(:)];

    rng(42);

    rSites = surfaceSufCond3D(rSites, F.c.CC, F.c.R);
    rSites = surfaceSufCond3D(rSites, F.f.pts, F.f.Gs);

    %% Build grid
    sites = [F.f.pts; rSites];

    eps_clip = 1e-8;

    sites(:,1) = min(max(sites(:,1), ...
        min(bdr(:,1))+eps_clip), max(bdr(:,1))-eps_clip);

    sites(:,2) = min(max(sites(:,2), ...
        min(bdr(:,2))+eps_clip), max(bdr(:,2))-eps_clip);

    sites(:,3) = min(max(sites(:,3), ...
        min(bdr(:,3))+eps_clip), max(bdr(:,3))-eps_clip);

    G = mirroredPebi3D(sites,bdr);
    G = computeGeometry(G);

    %% Compute robust geometric cell centroids
    geom_cc = zeros(G.cells.num,3);

    for ci = 1:G.cells.num
        [geom_cc(ci,:), ~] = computeRobustCellCentroid(G,ci);
    end

    G.cells.centroids = geom_cc;

    %% Classify cells by region
    c = G.cells.centroids;

    cl = c(:,2) > ...
        faultHeight2z(c(:,[1,3])) ...
        + faultHeight2x(c(:,[1,3])) - c(:,3);

    cr = c(:,2) < ...
        faultHeight1z(c(:,[1,3])) ...
        + faultHeight1x(c(:,[1,3])) - c(:,3);

    cm = ~cl & ~cr;

    %% Identify fault faces and adjacent cells
    region = zeros(G.cells.num,1);

    region(cl) = 1;
    region(cm) = 2;
    region(cr) = 3;

    neighbors = G.faces.neighbors;

    internalFaces = find(all(neighbors > 0,2));

    c1 = neighbors(internalFaces,1);
    c2 = neighbors(internalFaces,2);

    region1 = region(c1);
    region2 = region(c2);

    regionPairs = sort([region1,region2],2);

    % Fault 1: interface between regions 1 and 2
    fault1Mask = all(regionPairs == [1 2],2);

    % Fault 2: interface between regions 2 and 3
    fault2Mask = all(regionPairs == [2 3],2);

    fault1Faces = internalFaces(fault1Mask);
    fault2Faces = internalFaces(fault2Mask);

    fault1Cells = unique(neighbors(fault1Faces,:));
    fault2Cells = unique(neighbors(fault2Faces,:));

    faultCells = unique([fault1Cells;fault2Cells]);

    fprintf('\n=== Fault-adjacent cells ===\n');
    fprintf('Total cells: %d\n',G.cells.num);
    fprintf('Fault 1 faces: %d\n',numel(fault1Faces));
    fprintf('Fault 2 faces: %d\n',numel(fault2Faces));
    fprintf('Fault-adjacent cells: %d\n',numel(faultCells));

    %% Assign fault-adjacent cell indicator
    
    faultAdjacent = zeros(G.cells.num, 1);
    faultAdjacent(faultCells) = 1;
    
    % Store information in MRST grid
    G.cells.faultAdjacent = faultAdjacent;
    
    % Optional: store fault face IDs for later contact assembly
    G.faces.fault1 = false(G.faces.num, 1);
    G.faces.fault2 = false(G.faces.num, 1);
    
    G.faces.fault1(fault1Faces) = true;
    G.faces.fault2(fault2Faces) = true;
    
    %% Visualization: mesh only
    
    figure('Color', 'w');
    
    plotGrid(G, ...
        'FaceColor', [0.85 0.85 0.85], ...
        'FaceAlpha', 0.15, ...
        'EdgeColor', [0.3 0.3 0.3], ...
        'EdgeAlpha', 0.1);
    
    axis equal;
    view(120, 30);
    
    xlabel('x');
    ylabel('y');
    zlabel('z');
    
    title('Two-Fault Mesh');
end


%% Helper functions

function [centroid,volume] = computeRobustCellCentroid(G,ci)

    faceIdx = G.cells.faces( ...
        G.cells.facePos(ci):G.cells.facePos(ci+1)-1,1);

    ref = mean(G.faces.centroids(faceIdx,:),1);

    v_total = 0;
    c_weighted = [0,0,0];

    for k = 1:numel(faceIdx)
        f = faceIdx(k);

        sgn = 2*(G.faces.neighbors(f,1)==ci)-1;

        nStart = G.faces.nodePos(f);
        nEnd = G.faces.nodePos(f+1)-1;

        fnodes = G.faces.nodes(nStart:nEnd);
        nv = numel(fnodes);

        fc = mean(G.nodes.coords(fnodes,:),1);

        for j = 1:nv
            a = G.nodes.coords(fnodes(j),:) - ref;
            b = G.nodes.coords(fnodes(mod(j,nv)+1),:) - ref;
            cc = fc - ref;

            if sgn < 0
                [a,b] = deal(b,a);
            end

            vi = dot(a,cross(b,cc))/6;

            v_total = v_total + vi;
            c_weighted = c_weighted + vi*(a+b+cc)/4;
        end
    end

    volume = v_total;
    centroid = ref + c_weighted/v_total;
end


function R = circumradiusPerVertex(pts,conn,gamma)

    nPts = size(pts,1);

    p1 = pts(conn(:,1),:);
    p2 = pts(conn(:,2),:);
    p3 = pts(conn(:,3),:);

    L1 = sqrt(sum((p2-p1).^2,2));
    L2 = sqrt(sum((p3-p2).^2,2));
    L3 = sqrt(sum((p1-p3).^2,2));

    s = (L1+L2+L3)/2;

    area = sqrt(max(s.*(s-L1).*(s-L2).*(s-L3),0));

    Rcirc = (L1.*L2.*L3)./max(4*area,eps);

    verts = [conn(:,1);conn(:,2);conn(:,3)];
    rvals = [Rcirc;Rcirc;Rcirc];

    R = accumarray(verts,rvals,[nPts,1],@max);

    R = R/gamma;
end