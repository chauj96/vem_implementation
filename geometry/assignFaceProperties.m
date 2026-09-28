function face_struct = assignFaceProperties(face_struct)
% identify boundary faces and tag boundary conditions
%
% manufactured solution for spatial convergence test
%   Mechanics : prescribed displacement on entire boundary -> u = u_exact
%   Flow      : prescribed pressure on entire boundary -> p = p_exact

    nFaces = numel(face_struct);

    for f = 1:nFaces

        % boundary flag
        face_struct(f).is_boundary = false;

        % mechanics
        face_struct(f).is_displacement_bc = false;
        face_struct(f).is_traction_bc = false;

        % flow
        face_struct(f).is_pressure_bc = false;
        face_struct(f).is_flux_bc = false;

        % exterior face
        if isscalar(face_struct(f).cells)

            face_struct(f).is_boundary = true;

            % Gamma_u = boundary of Omega
            face_struct(f).is_displacement_bc = true;

            % Gamma_p = boundary of Omega
            face_struct(f).is_pressure_bc = true;

        end
    end

end
