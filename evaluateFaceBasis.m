function phi = evaluateFaceBasis(j,x,xf,Qf,t1,t2,n,Af,m20,m02,m11)

    local = Qf' * (x - xf);

    xt = local(1);
    yt = local(2);

    switch j

        case 1
            phi = t1/Af;

        case 2
            phi = t2/Af;

        case 3
            phi = n/Af;

        case 4
            phi = (-yt*t1 + xt*t2)/(Af*sqrt(m20+m02));

        case 5
            phi = (xt*n)/(Af*sqrt(m20));

        case 6
            phi = (yt*n)/(Af*sqrt(m02));

    end

end