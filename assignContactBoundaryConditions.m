function face_struct = assignContactBoundaryConditions(face_struct, sigma0)
% define component-wise boundary conditions for the inclined fracture compression problem, including uniaxial compression,
% traction-free lateral boundaries, and plane-strain constraints.
%
% store displacement and traction BC metadata on boundary faces
% without imposing constraints on the global system.
%
% sigma0: positive compression magnitude (MPa)

    assert(sigma0 > 0, 'Compression magnitude must be positive.');

    for f = 1:numel(face_struct)

        % initialize component-wise BC information
        face_struct(f).displacement_mask = false(3,1);
        face_struct(f).displacement_value = zeros(3,1);

        face_struct(f).traction_mask = false(3,1);
        face_struct(f).traction_value = zeros(3,1);

        if ~face_struct(f).is_boundary
            continue;
        end

        side = face_struct(f).boundary_side;

        switch side
        
            case 'xmin'
                % uniaxial compression
                face_struct(f).traction_mask(:) = true;
                face_struct(f).traction_value = [sigma0; 0; 0];
        
            case 'xmax'
                % uniaxial compression
                face_struct(f).traction_mask(:) = true;
                face_struct(f).traction_value = [-sigma0; 0; 0];
        
            case {'ymin','ymax'}
                % traction-free boundaries
                face_struct(f).traction_mask(:) = true;
        
            case {'zmin','zmax'}
                % plane strain
                face_struct(f).displacement_mask(3) = true;
                face_struct(f).traction_mask([1,2]) = true;
        
            otherwise
                error('Unknown boundary side on face %d.', f);
        
        end

        face_struct(f).is_displacement_bc = any(face_struct(f).displacement_mask);

        face_struct(f).is_traction_bc = any(face_struct(f).traction_mask);

        assert(all(face_struct(f).displacement_mask | face_struct(f).traction_mask), 'Incomplete boundary conditions on face %d.', f);

        assert(~any(face_struct(f).displacement_mask & face_struct(f).traction_mask), 'Conflicting boundary conditions on face %d.', f);

    end

    fprintf('\n[# CONTACT BOUNDARY CONDITIONS]\n');
    fprintf('  Compression magnitude: %.2f MPa\n', sigma0);
    fprintf('  Boundary faces: %d\n', nnz([face_struct.is_boundary]));

end
