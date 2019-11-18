%This function accepts the fin thickness, the minor diameter of the
%elliptic fin, the outer radius of the cylinder, the thermal conductivity
%of the fins, and the air-side HTC, and returns the optimum major radius of
%the ellipse that will return the highest fin effectiveness. 

function r1 = optimumFin(t,r2,rb,k,h)
    Bi = h*rb/k;
    r1 = 0;
    r12 = r2+0.05;
    conv = abs(r1-r12);
    while conv > 0.001
        r1 = r12;
        fr = func(r1);
        dfr = dfunc(r1);
        r12 = r1-fr/dfr;
        conv = abs(r1-r12);
    end

    function fr = func(r1)
        zeta=t/rb;
        R1=r1/rb;
        R2=r2/rb;
        U=(R1*R2-1)*zeta;
        fr=0.447*(U*Bi)^0.317*(U)^(0.22*R2)-zeta;
    end

    function dfr = dfunc(r1)
        U=t*(r1*r2-rb^2)/rb^3;
        dU=t*r2/rb^3;
        R=r2/rb;
        term1=0.447*0.317*(dU*Bi)*(U*Bi)^-0.683*U^(0.22*R);
        term2=0.447*0.22*R*dU*(U*Bi)^0.317*U^(0.22*R-1);
        dfr=term1+term2;
    end
end