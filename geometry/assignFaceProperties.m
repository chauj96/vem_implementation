function face_struct = assignFaceProperties(face_struct)
% identify boundary faces and tag boundary conditions
%
% Consolidation of a soil layer
%
% z = 0 : top
%   Mechanics : prescribed traction
%   Flow      : prescribed pressure (p = 0)
%
% z = 1 : bottom
%   Mechanics : prescribed normal displacement (u_z = 0)
%   Flow      : prescribed normal flux (w.n = 0)
%
% Side boundaries
%   Mechanics : prescribed normal displacement (u.n = 0)
%   Flow      : prescribed normal flux (w.n = 0)

    nFaces = numel(face_struct);

    h = 1.0;
    tol = 1e-10;

    for f = 1:nFaces

        % boundary flag
        face_struct(f).is_boundary = false;

        % boundary location
        face_struct(f).is_top = false;
        face_struct(f).is_bottom = false;
        face_struct(f).is_side = false;

        % mechanics
        face_struct(f).is_displacement_bc = false;
        face_struct(f).is_traction_bc = false;

        % flow
        face_struct(f).is_pressure_bc = false;
        face_struct(f).is_flux_bc = false;

        % exterior face
        if isscalar(face_struct(f).cells)

            face_struct(f).is_boundary = true;

            zf = face_struct(f).center(3);

            if abs(zf) < tol

                % Top boundary (z = 0)
                face_struct(f).is_top = true;

                face_struct(f).is_traction_bc = true;
                face_struct(f).is_pressure_bc = true;

            elseif abs(zf - h) < tol

                % Bottom boundary (z = h)
                face_struct(f).is_bottom = true;

                face_struct(f).is_displacement_bc = true;
                face_struct(f).is_flux_bc = true;

            else

                % Side boundaries
                face_struct(f).is_side = true;

                face_struct(f).is_displacement_bc = true;
                face_struct(f).is_flux_bc = true;

            end
        end
    end
end
