function face_struct = assignFaceProperties(face_struct)
% identify exterior boundary faces and assign boundary condition tags.
% we have linear elasticity with an internal contact interface.
%
% patch test:
%   prescribed displacement on the entire exterior boundary.
%
%   Gamma_D = boundary(Omega)
%   Gamma_N = empty
%
% internal fracture faces are not exterior boundary faces.
% they are handled separately by the contact formulation.

    nFaces = numel(face_struct);

    for f = 1:nFaces

        % initialize boundary flags
        face_struct(f).is_boundary = false;

        % mechanical boundary conditions
        face_struct(f).is_displacement_bc = false;
        face_struct(f).is_traction_bc = false;

        % exterior boundary faces have one adjacent cell
        if isscalar(face_struct(f).cells)

            face_struct(f).is_boundary = true;

            % prescribed displacement on all exterior faces
            face_struct(f).is_displacement_bc = true;

        end

    end

    fprintf('\n[# BOUNDARY FACE PROPERTIES]\n');

    nBoundary = sum([face_struct.is_boundary]);
    nDirichlet = sum([face_struct.is_displacement_bc]);
    nNeumann = sum([face_struct.is_traction_bc]);

    fprintf('  Exterior boundary faces: %d\n', nBoundary);
    fprintf('  Displacement BC faces:   %d\n', nDirichlet);
    fprintf('  Traction BC faces:       %d\n', nNeumann);

end

