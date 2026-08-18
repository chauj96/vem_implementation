function face_struct = assignFacePropertiesDirichlet(face_struct)
% tag every exterior face as Dirichlet
%
% variant of assignFaceProperties for the convergence tests, where boundary
% conditions are taken from the analytical solution on the whole boundary.
% in the dual mixed formulation a prescribed displacement is a natural
% condition, so no stress DOF is eliminated and assembleNeumann is not needed

    nFaces = numel(face_struct);

    for f = 1:nFaces

        onBoundary = isscalar(face_struct(f).cells);

        face_struct(f).is_boundary = onBoundary;
        face_struct(f).is_dirichlet = onBoundary;
        face_struct(f).is_neumann = false;

    end

end
