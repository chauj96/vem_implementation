function cellData = computeCellFields(cell_struct, R_local, sigma_h, u_h, errors)
% build the cell-centred fields written to the VTU file
%
% u_h|_E = alpha + omega x (x - x_E), so alpha = u_h(x_E) and omega is the
% infinitesimal rotation. the stress tensor is S(c) with c = R_E sigma_h|_E

    nCells = numel(cell_struct);

    displacement = zeros(nCells,3);
    rotation = zeros(nCells,3);
    stress = zeros(nCells,9);
    vonMises = zeros(nCells,1);

    % I1 = tr(S),  I2 = (I1^2 - tr(S^2))/2,  I3 = det(S)
    I1 = zeros(nCells,1);
    I2 = zeros(nCells,1);
    I3 = zeros(nCells,1);

    s = 1/sqrt(2);

    for e = 1:nCells

        idx = 6*(e-1) + (1:6);

        displacement(e,:) = u_h(idx(1:3)).';
        rotation(e,:) = u_h(idx(4:6)).';

        % 6(f-1)+i for the faces of E
        faces = cell_struct(e).faces(:).';
        stressDofs = reshape(6*(faces - 1) + (1:6).', 1, []);

        c = R_local{e} * sigma_h(stressDofs);

        S = [   c(1), s*c(6), s*c(5); ...
             s*c(6),    c(2), s*c(4); ...
             s*c(5), s*c(4),    c(3)];

        % S = S', so column-major and row-major agree
        stress(e,:) = S(:).';

        trS = trace(S);

        % S = S'  =>  tr(S^2) = sum_ij S_ij^2
        I1(e) = trS;
        I2(e) = 0.5 * (trS^2 - sum(S(:).^2));
        I3(e) = det(S);

        % sqrt(3 J2),  J2 = |dev(S)|^2 / 2
        dev = S - (trS/3)*eye(3);

        vonMises(e) = sqrt(1.5 * sum(dev(:).^2));

    end

    cellData.displacement = displacement;
    cellData.rotation = rotation;
    cellData.stress = stress;
    cellData.stress_I1 = I1;
    cellData.stress_I2 = I2;
    cellData.stress_I3 = I3;
    cellData.von_mises = vonMises;

    % |E| and h_E, useful for spotting degenerate cells
    cellData.volume = [cell_struct.volume].';
    cellData.diameter = [cell_struct.diameter].';

    if nargin >= 5 && isstruct(errors) && isfield(errors, 'u_exact')

        % |alpha_h - alpha| per cell
        D = reshape(u_h - errors.u_exact, 6, nCells);

        cellData.displacement_error = vecnorm(D(1:3,:), 2, 1).';

    end

end
