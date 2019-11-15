%This function accepts the fin thickness, the minor diameter of the
%elliptic fin, the outer radius of the cylinder, the thermal conductivity
%of the fins, and the air-side HTC, and returns the optimum major radius of
%the ellipse that will return the highest fin effectiveness. 

function ra = optimumFin(t,rb,rc,k,h)
    Bi = h*rc/k;
    ra = 0;
    ra2 = rb+0.5;
    conv = abs(ra-ra2);
    while conv > 0.001
        ra = ra2;
        fr = func(ra);
        dfr = dfunc(ra);
        ra2 = ra-fr/dfr;
        conv = abs(ra-ra2);
    end

    function fr = func(ra)
        U=t*(ra*rb-rc^2)/rc^3;
        fr=0.447*(U*Bi)^0.317*(U)^(0.22*rb)-t/rb;
    end

    function dfr = dfunc(ra)
        U=t*(ra*rb-rc^2)/rc^3;
        dU=t*rb/rc^3;
        term1=0.447*0.317*(dU*Bi)*(U*Bi)^-0.683*U^(0.22*rb);
        term2=0.447*0.22*rb*dU*(U*Bi)^0.317*U^(0.22*rb-1);
        dfr=term1+term2;
    end
end