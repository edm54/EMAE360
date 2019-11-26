%{
This program finds the required velocity for air-cooling cylinders. From
this, the required power to accelerate this amount of air is calculated as
well. 
%}
clear
eq=1;

b = 0.0705; %m, bore
s = 0.0641; %m, stroke
Ag=2*pi*s*b;
thick = 0.017381; %m, cylinder thickness, confirm with team

Tw=600; %K, a set value. 
To=300; %K, atmospheric temperature

%Fin characteristics 
tf=0.0025; %Thickness of fins
t=tf/2; %Half thickness of fins, used to calculate m
S=0.003; %fin spacing
k=151;

%Elliptic Fin Characteristics
term = 0.0707/2;
rc = b/2; %Inside of cylinder
r2=rc+thick; %Outside Radius of Cylinder 
rb = r2+term; %minor axis radius of fins
N1=floor((0.12)/(S+tf))+1; %stroke plus thickness of cylinder

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
    hflux(i)=AverageTempAndH(nv(i),eq,rpm,Tw,0)*10^6; %W/m^2
    Q=hflux(i)*Ag;
    Rw=log(r2/rc)/(2*pi*s*k);
    Ts1=Tw-Q*Rw;

    v = 3; %m/s
    Ts2 = 1000; %arbitrary
    while Ts2>Ts1
        v = v+1;
        hc = CoolingHTC(K,B,v,t,S); %Heat Transfer Coefficient between Fin and air
        k = 151; %Thermal conductivity of fins, W/(m K)
        m = sqrt(hc/(k*t));
        %ra = optimumFin(t,rb,r2,k,hc);
        ra =0.14;
        nfe = FinEffEll(ra,rb,r2,m);
        Afe=EllipticArea(ra,rb,r2,N1);
%        Afe=2*N1*pi*(ra^2-r2^2);
        Ate=Afe+CylinArea(r2,S,N1);
        Rfe=1/(nfe*hc*Afe);
        Rte=1/(hc*(Ate-Afe));
        Re=Rte*Rfe/(Rte+Rfe);
        Ts2=Q*Re+To;
    end
    vreq(i) = v*2.237;
    %fprintf("Required airspeed is %d mph. \n", v)
end
plot(N,vreq)
title("Required Air Velocity for Cooling vs RPM")
xlabel("RPM")
ylabel("Air Velocity, mph")
grid on

function nf = FinEffEll(ra,rb,rc,m)

    %{
    %Bessel function variant for circular annular fins
    term1 = 2*rc / (m*(ra^2-rc^2));
    term2 = besselk(1,m*rc)*besseli(1,m*ra)-besseli(1,m*rc)*besselk(1,m*ra);
    term3 = besselk(0,m*rc)*besseli(1,m*ra)+besseli(0,m*rc)*besselk(1,m*ra);
    nf = term1*(term2/term3);
    %}

    %Elliptic Annular fin efficiency
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