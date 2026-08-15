function D_local = computeDE(cell_struct, A_local, B_local)
% construct the local divergence reconstruction matrix D_E
%
%   D_E = [ B_E(1:3,:) / |E| ;  A_E \ B_E(4:6,:) ]
%
% one factorization of A_E for all columns

    nCells = numel(A_local);
    D_local = cell(nCells, 1);

    for e = 1:nCells

        A_E = A_local{e};
        % assert(rcond(A_E) > 1e-12, 'A_E is nearly singular.');

        volume = cell_struct(e).volume;
        B_E = B_local{e}.matrix;

        % translation
        alpha = B_E(1:3, :) / volume;

        % rotation
        omega = A_E \ B_E(4:6, :);

        D_local{e} = [alpha; omega];

    end

end
