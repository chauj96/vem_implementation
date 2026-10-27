function K = assembleContactGlobalSystem(A, B, Ac, T, J)
% assemble the fully mixed elasticity-contact system.
%
% unknown ordering:
%   x = [sigma; u; tc; uc]
%
% block system:
%   [ A   B'   0   T' ]
%   [ B    0   0    0 ]
%   [ 0    0   Ac   J ]
%   [ T    0   J'   0 ]

    nSigma = size(A, 1);
    nU = size(B, 1);
    nTc = size(Ac, 1);
    nUc = size(T, 1);

    % Check operator dimensions
    assert(isequal(size(A), [nSigma, nSigma]));
    assert(isequal(size(B), [nU, nSigma]));
    assert(isequal(size(Ac), [nTc, nTc]));
    assert(isequal(size(T), [nUc, nSigma]));
    assert(isequal(size(J), [nTc, nUc]));

    % Sparse zero blocks
    Zuu = sparse(nU, nU);
    Zst = sparse(nSigma, nTc);
    Zut = sparse(nU, nTc);
    Zuc = sparse(nU, nUc);
    Zts = sparse(nTc, nSigma);
    Ztu = sparse(nTc, nU);
    Zus = sparse(nUc, nU);
    Zcc = sparse(nUc, nUc);

    % Assemble global block matrix
    K = [ A,    B',   Zst,  T';
          B,    Zuu,  Zut,  Zuc;
          Zts,  Ztu,  Ac,   J;
          T,    Zus,  J',   Zcc ];

end
