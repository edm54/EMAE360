%{
This program finds the required velocity for aircooling cylinders. From
this, the required power to accelerate this amount of air is calculated as
well. 
%}
clear
eq=1;

b = 0.0705; %m, bore
s = 0.0641; %m, stroke
Ag=2*pi*s*b;
thick = 0.017381; %m, cylinder thickness, confirm with team

Tw=530; %K, a set value. 
To=300; %K, atmospheric temperature

%Fin characteristics 
tf=0.0025; %Thickness of fins
t=tf/2; %Half thickness of fins, used to calculate m
S=0.003; %fin spacing
k=151;

%Elliptic Fin Characteristics
term = (0.0707-0.001)/2;
rc = b/2; %Inside of cylinder
r2=rc+thick; %Outside Radius of Cylinder 
rb = 0.0705/2+term; %minor axis radius of fins
N1=floor((s+0.0487)/(S+tf))+1; %stroke plus thickness of cylinder

%Rectangular Fin Characteristics
x = 0.06;
L2 = 0.04922;
L1 = 0.042;
N2 = 2*floor(L1/(S+tf))+floor(L2/(S+tf)); %includes both sides of the cylinder

if S/b >= 0.14
    K = 0.36;
    B = 0.55;
else
    K = 0.62;
    B = 0.27;
end
N=800:100:9300;
nv = volumetric_efficiency(N);
close all
vreq = zeros(length(N),1);
hflux = zeros(length(N),1);
for i = 1:length(N)
    rpm = N(i);
    hflux(i)=AverageTempAndH(nv(i),eq,rpm,0)*10^6; %W/m^2
    Q=hflux(i)*Ag;
    Rw=log(r2/rc)/(2*pi*s*k);
    Ts1=Tw-Q*Rw;

    v = 3; %m/s
    Ts2 = 1000; %arbitrary
    while Ts2>Ts1
        v = v+0.01;
        hc = CoolingHTC(K,B,v,t,S); %Heat Transfer Coefficient between Fin and air
        k = 151; %Thermal conductivity of fins, W/(m K)
        m = sqrt(hc/(k*t));
        ra = optimumFin(t,rb,r2,k,hc);
        nfe = FinEffEll(ra,rb,r2,m);
        nfr = FinEffRect(x,t,m);
        Afe=EllipticArea(ra,rb,r2,N1);
        Afr=RectArea(x,N2);
        Ate=Afe+CylinArea(r2,S,N1);
        Atr=Afr+SurfArea(S,N2);
        Rfe=1/(nfe*hc*Afe);
        Rfr=1/(nfr*hc*Afr);
        Rte=1/(hc*(Ate-Afe));
        Rtr=1/(hc*(Atr-Afr));
        Rf=Rfe*Rfr/(Rfe+Rfr);
        Re=Rte*Rfe/(Rte+Rfe);
        Rr=Rtr*Rfr/(Rtr+Rfr);
        R1=Re*Rr/(Re+Rr);
        Ts2=Q*R1+To;
    end
    vreq(i) = v*2.237;
    %fprintf("Required airspeed is %d mph. \n", v)
end
plot(N,vreq)
title("Required Air Velocity for Cooling vs RPM")
xlabel("RPM")
ylabel("Air Velocity, mph")

function nf = FinEffEll(ra,rb,rc,m)
%{
    Bessel function variant for circular annular fins
    term1 = 2*r1 / (m*(r2^2-r1^2));
    term2 = besselk(1,m*r1)*besseli(1,m*r2)-besseli(1,m*r1)*besselk(1,m*r2);
    term3 = besselk(0,m*r1)*besseli(1,m*r2)+besseli(0,m*r1)*besselk(1,m*r2);
    nf = term1*(term2/term3);
%}
    L=((ra-rc)+(rb-rc))/2;
    Rf=sqrt(ra*rb/rc^2);
    psi=1+0.17912*log(Rf);
    nf=(tanh(m*psi*L)/(m*psi*L))^psi;
end

function nf = FinEffRect(x,t,m)
    Lc = x+t;
    nf = tanh(m*Lc)/(m*Lc);
end

function hc = CoolingHTC(K,B,v,tf,S)
    a=0.07;
    k=.02435; %thermal conductivity of air, W/(m K)
    b=0.0705; %bore, m
    mu=18.5e-06; %Pa*s
    Re=v*b/mu;
    X=(tf/S+1)^0.55 * (1-K*Re^(-a)/((S/b)^B))^0.55;
    Nu=0.446*X*Re^0.55;
    hc=Nu*k/b;
end

function Afe = EllipticArea(ra,rb,rc,N1)
    Afe = 2*pi*(ra*rb-rc^2)*N1;
end

function Afr = RectArea(x,N2)
    L =0.1412; %m
    Afr = 2*L*x*N2;
end

function Afc = CylinArea(rc,S,N)
    Afc = 2*pi*rc*S*(N-1);
end

function Atr = SurfArea(S,N)
    Atr = 2*S*0.1412*N;
end