function A_local = computeAE_Tonon(cell_struct, face_struct, V3)
% Construct inertia matrix A_E analytically by tetrahedral decomposition,
% referred to F.Tonon (2005)

    nCells = numel(cell_struct);
    A_local = cell(nCells,1);

    for e = 1:nCells

        xE = cell_struct(e).center(:);
        A_E = zeros(3,3);
        cell_faces = cell_struct(e).faces;

        for lf = 1:numel(cell_faces)

            f = cell_faces(lf);
            xf = face_struct(f).center(:);
            verts = face_struct(f).verts;
            nVerts = numel(verts);

            for k = 1:nVerts

                kp1 = mod(k, nVerts)+1;

                vi = V3(verts(k), :)';
                vip1 = V3(verts(kp1), :)';

                % relative vectors
                a = xf - xE;
                b = vi - xE;
                c = vip1 - xE;

                % Jacobian
                B = [a b c];
                J = abs(det(B));

                % Tonon's closed-form inertia matrix
                S = a*a' + b*b' + c*c' + (a+b+c) * (a+b+c)';
                A_T = (J/120) * (trace(S) * eye(3) - S);
                A_E = A_E + A_T;
            end
        end

        A_local{e} = A_E;

    end

end