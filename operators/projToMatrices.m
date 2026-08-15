function matrices = projToMatrices(P_local)
% convert the projection coefficients into symmetric stress tensors

    numBasisFunctions = size(P_local,2);
    matrices = cell(numBasisFunctions,1);
    
    for j = 1:numBasisFunctions
      matrices{j} = ...
        [P_local(1,j),P_local(6,j)/sqrt(2),P_local(5,j)/sqrt(2);
        P_local(6,j)/sqrt(2),P_local(2,j),P_local(4,j)/sqrt(2);
        P_local(5,j)/sqrt(2),P_local(4,j)/sqrt(2),P_local(3,j)];
    end

end
