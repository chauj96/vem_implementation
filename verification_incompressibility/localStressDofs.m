function sl = localStressDofs(e, cell_struct, face_struct, face_global_geom, B_local, sigma_h)
% sigma_h restricted to element e, converted to the element's own face
% convention: x_local = T_E x_global

    fl = cell_struct(e).faces(:).';
    sl = zeros(6*numel(fl), 1);

    for lf = 1:numel(fl)
        f = fl(lf);
        d = faceSignTransformation(f, B_local{e}.geom(lf), face_global_geom(f));
        sl(6*(lf-1) + (1:6)) = d .* sigma_h(6*(f-1) + (1:6));
    end

end
