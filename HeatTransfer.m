%{
This program finds the required velocity for air-cooling cylinders. From
this, the required power to accelerate this amount of air is calculated as
well. 
%}
clear
eq=1;
N=800:100:9300;
N=9300;
Wt=zeros(length(N),1);
Power=zeros(length(N),1);
P_hp=zeros(length(N),1);
vel=zeros(length(N),1);
vreq = zeros(length(N),1);
hflux = zeros(length(N),1);
Ts2 = zeros(length(N),1);
Tw = zeros(length(N),1);

af = 14.7;
Vtdc = 2.781e-5; %m^3
V1 = Vtdc*10;
Qlhv = 44e6; %j/kg from Heywood App D
p1 = 101325;
RPM = [2100 15000];
mech_eff = [ .9 .65]; 
[nv,pump_work] = volumetric_efficiency(N);
for i = 1:length(N)
    ma = 3.335e-4; %kg
    close all
    ma = ma*nv(i);
    mf = ma/af;
    Qin = mf*combustionEff(1)*Qlhv/(p1*V1); 
    net_work = 6 * (FiniteHeatWoschni(Qin, N(i), ma, 0)); %kJ   
    if N(i)<= 2100
        mechanical_eff = .9;
    else
        mechanical_eff = interp1(RPM, mech_eff, N(i), 'linear');
    end
    efficiency = mechanical_eff;
    Wt(i) = efficiency * net_work/6  - pump_work(i)/1000;

    Power(i) = 6 * Wt(i) * N(i)/120; % watts
    P_hp(i) = Power(i)/0.745699872;
    vel(i) = velocity_finder(P_hp(i)); %mph
end

b = 0.0705; %m, bore
s = 0.0641; %m, stroke
Ag=2*pi*s*b;
thick = 0.017381; %m, cylinder thickness, confirm with team

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
close all
v = 3; %m/s arbitrary
for i = 1:length(N)
    rpm = N(i);
    if i == 1
        Tw(i)=530; %K, a set value. 
    else
        Tw(i)=Tw(i-1);
    end
    hflux(i)=AverageTempAndH(nv(i),eq,rpm,Tw(i),0)*10^6; %W/m^2
    Q=hflux(i)*Ag;
    Rw=log(r2/rc)/(2*pi*s*k);
    Ts1=Tw(i)-Q*Rw;
    Ts2(i) = 1000; %arbitrary
    while Ts2(i)>Ts1   
        if (2.237*v)>=vel(i)
            Tw(i) = Tw(i)+1;
            hflux(i)=AverageTempAndH(nv(i),eq,rpm,Tw(i),0)*10^6; %W/m^2
            Q=hflux(i)*Ag;
            Rw=log(r2/rc)/(2*pi*s*k);
            Ts1=Tw(i)-Q*Rw;
        else
            v = v+0.05;
        end
        hc = CoolingHTC(K,B,v,t,S); %Heat Transfer Coefficient between Fin and air
        k = 151; %Thermal conductivity of fins, W/(m K)
        m = sqrt(hc/(k*t));
        %ra = optimumFin(t,rb,r2,k,hc);
        ra =0.1201;
        nfe = FinEffEll(ra,rb,r2,m);
        Afe=EllipticArea(ra,rb,r2,N1);
%        Afe=2*N1*pi*(ra^2-r2^2);
        Ate=Afe+CylinArea(r2,S,N1);
        Rfe=1/(nfe*hc*Afe);
        Rte=1/(hc*(Ate-Afe));
        Re=Rte*Rfe/(Rte+Rfe);
        Ts2(i)=Q*Re+To;
    end
    vreq(i) = v*2.237;
    %fprintf("Required airspeed is %d mph. \n", v)
end
hold on
figure
plot(N,vreq,N,vel)
title("Required Air Velocity for Cooling vs RPM")
xlabel("RPM")
ylabel("Air Velocity, mph")
legend("Required Air Velocity", "Max motorcycle Veloity")
grid on

figure
plot(N,Tw,N,Ts2)
ylabel("Temperature, K")
xlabel("RPM")
legend("Internal Wall temperature", "Cylinder Surface Temperature")
hold off

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