function [neuDofs, neuVals] = computeNeumannData(face_struct, face_global_geom, lambda, mu)
% compute prescribed Neumann stress DOFs in the global/reference face convention

    % exact stress corresponding to u = [x+1; 0; 0]
    sigma = [2*mu + lambda, 0, 0;
             0, lambda, 0;
             0, 0, lambda];

    neumannFaces = find([face_struct.is_neumann]);

    nNeumannFaces = numel(neumannFaces);
    nNeumannDofs = 6 * nNeumannFaces;

    neuDofs = zeros(1,nNeumannDofs);
    neuVals = zeros(nNeumannDofs,1);

    for k = 1:nNeumannFaces

        f = neumannFaces(k);

        Af = face_struct(f).area;

        geom = face_global_geom(f);

        ng = geom.n;
        t1g = geom.t1;
        t2g = geom.t2;

        % prescribed traction in the global face convention
        traction = sigma * ng;

        % coefficients in the global face traction basis
        coeff = [ ...
            Af * dot(traction,t1g);
            Af * dot(traction,t2g);
            Af * dot(traction,ng);
            0;
            0;
            0];

        globalIdx = 6*(f-1) + (1:6);
        localIdx = 6*(k-1) + (1:6);

        neuDofs(localIdx) = globalIdx;
        neuVals(localIdx) = coeff;

    end

end