%{
This program finds the required velocity for aircooling cylinders. From
this, the required power to accelerate this amount of air is calculated as
well. 
%}
clear
b = 0.0705; %m, bore
s = 0.0641; %m, stroke
thick = 0.01; %m, cylinder thickness, confirm with team

Tw=530; %K, a set value. 
To=300; %K, atmospheric temperature
hflux=0.3165e6; %W/m^2

%Fin characteristics 
L=0.12/2; %cylinder spacing divided by two
tf=0.006; %Thickness of fins
t=tf/2; %Half thickness of fins, used to calculate m
S=0.006; %fin spacing
k=151;
N=floor(b/(S+tf))+1;

term = (0.0707-0.002)/2;
rc = b/2;
rb = 0.0705/2+term;
ra = optimumFin(t,rb,b/2,k,hc);

r1=b/2 + thick; %Outside of cylinder
r2=r1+L; %Outside Radius of Fins 
Ag=2*pi*s*b;
Af=(2*pi*(r2^2-r1^2)+2*pi*r2*t)*2;
At=(N*Af+2*pi*r1*S*(N-1));

Q=hflux*Ag;
Rw=log(r2/r1)/(2*pi*s*k);

Ts1=Tw-Q*Rw;
Ts2=0;
v = 6; %m/s

if S/b >= 0.14
    K = 0.36;
    B = 0.55;
else
    K = 0.62;
    B = 0.27;
end
while abs(Ts1-Ts2) > 0.001
    v = v+0.0001;
    hc = CoolingHTC(K,B,v,t,S); %Heat Transfer Coefficient between Fin and air
    k = 151; %Thermal conductivity of fins, W/(m K)
    m = sqrt(2*hc/(k*t));
    ra = optimumFin(t,rb,rc,k,hc);
    nf = FinEff(ra,rb,rc,m);
    Af=2*pi*(ra*rb-rc^2)+2*pi*r2*tf;
    At=N*Af+2*pi*r1*S*(N-1);
    Rf=1/(N*nf*hc*Af);
    Rt=1/(hc*(At-N*Af));
    R1=Rf*Rt/(Rf+Rt);
    Ts2=Q*R1+To;
end

function nf = FinEff(ra,rb,rc,m)
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