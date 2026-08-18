function face_struct = assignFacePropertiesSlab(face_struct, tol)
% boundary tagging for the quasi-2D slab
%
% the slab's z-normal faces are an artefact of embedding a plane problem in 3D,
% not physical boundaries. leaving their stress DOFs free lets sigma_zz absorb
% lambda * tr(eps_h), which is O(lambda h^2) and does not decay until
% h^2 < 1/lambda. they are therefore given the exact traction as an essential
% (Neumann) condition, while the genuine x- and y-boundaries stay Dirichlet

    if nargin < 2 || isempty(tol), tol = 1e-10; end

    for f = 1:numel(face_struct)

        onBoundary = isscalar(face_struct(f).cells);

        n = face_struct(f).normal(:);
        isZ = abs(abs(n(3)) - 1) < tol;

        face_struct(f).is_boundary = onBoundary;
        face_struct(f).is_neumann = onBoundary && isZ;
        face_struct(f).is_dirichlet = onBoundary && ~isZ;

    end

end
