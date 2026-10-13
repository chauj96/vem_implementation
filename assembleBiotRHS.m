function [f_sigma, f_u, f_w, f_p] = assembleBiotRHS(cell_struct, face_struct, face_global_geom, t)

    nCells = numel(cell_struct);
    nFaces = numel(face_struct);

    f_sigma = zeros(6*nFaces, 1);
    f_u = zeros(6*nCells, 1);
    f_w = zeros(nFaces, 1);
    f_p = zeros(nCells, 1);

end

