function R = evaluateRMBasis(i, x, xE)

    r = x - xE;

    switch i

        case 1
            R = [1;0;0];

        case 2
            R = [0;1;0];

        case 3
            R = [0;0;1];

        case 4
            R = cross([1;0;0], r);

        case 5
            R = cross([0;1;0], r);

        case 6
            R = cross([0;0;1], r);

    end

end