function [t1, t2, Qf] = constructFaceFrame(X, n)
% construct an orthonormal local frame {t1,t2,n}.

    nVertices = size(X, 1);
    t1 = []; 

    for k = 1:nVertices

        kp1 = mod(k, nVertices) + 1;
        edge = X(kp1,:)' - X(k,:)';

        edge = edge - dot(edge, n) * n;

        if norm(edge) > 1e-14
            t1 = edge / norm(edge);
            break;
        end
    end

    if isempty(t1)
        error('Could not construct t1: all face edges are degenerate.');
    end

    t2 = cross(n, t1);
    t2 = t2 / norm(t2);

    % reorthogonalize t1 to reduce floating-point drift
    t1 = cross(t2, n);
    t1 = t1 / norm(t1);

    Qf = [t1, t2, n];
end