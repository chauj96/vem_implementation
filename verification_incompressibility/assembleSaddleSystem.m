function [K_global, B_global] = assembleSaddleSystem(cell_struct, face_struct, face_global_geom, B_local, K_cons, K_stab)
% K and B from assembleGlobalMatrices, independent of which output convention
% that function currently uses
%
%   3 outputs : [K_cons_global, K_stab_global, B_global]
%   4 outputs : [K_global, B_global, K_cons_global, K_stab_global]

    if nargout('assembleGlobalMatrices') >= 4
        [K_global, B_global] = assembleGlobalMatrices(cell_struct, face_struct, ...
            face_global_geom, B_local, K_cons, K_stab);
    else
        [Kc, Ks, B_global] = assembleGlobalMatrices(cell_struct, face_struct, ...
            face_global_geom, B_local, K_cons, K_stab);
        K_global = Kc + Ks;
    end

end
