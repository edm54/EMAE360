%Motored Engine
clear
r = 10;
T1 = 300; %K 
P1 = 101325; %Pa
V1 = 0.000275; %m^3
R = 287; %J/kg K

thetai = 0;
mstep = 0.001;
thetaf = 2*pi;
mprops.theta = thetai:mstep:thetaf;
mprops.theta = mprops.theta.';
mprops.vol = zeros(length(mprops.theta),1);
mprops.press = zeros(length(mprops.theta),1);
mprops.temp = zeros(length(mprops.theta),1);
mprops.rho = zeros(length(mprops.theta),1);
gamma = zeros(length(mprops.theta),1);

mprops.vol(1) = nondimV(thetai, r)*V1;
mprops.press(1) = 1;
mprops.temp(1) = T1;
k = P1*(V1)^(calc_gamma(T1));
mprops.rho(1) = P1/(287*T1);
ma = mprops.rho(1)*(mprops.vol(1));

for i = 2:length(mprops.theta)
    gamma(i) = calc_gamma(mprops.temp(i-1));
    [mprops.vol(i),~] = nondimV(mprops.theta(i),r);
    mprops.vol(i) = mprops.vol(i)*V1;
    mprops.rho(i) = ma / (mprops.vol(i));
    mprops.press(i) = k / (mprops.vol(i))^(gamma(i));
    mprops.temp(i) = mprops.press(i)/(mprops.rho(i)*R);
end

%Real Engine Cycle with Heat Transfer from the Cylinder


%%Overall heat Transfer coefficient

%% Real Engine Cycle with Heat Transfer from the Cylinde
pi = 3.142;
step = 0.001; %crankangle interval for calculation
thetai = 0; %initial crankangle, rad
thetaf =2*pi; %final crankangle througout steps
NN = (thetaf - thetai)/step + 1;

thetas = (8/9)*pi; %start of heat release (deg)
thetad = pi/3; %duration of heat release (deg)
num = thetas:step:(thetas + thetad);

r = 10; %compression ratio
ma = 0.000336; %kg, Get from intake/outtake program
af = 14.7; %Air to fuel ratio
mf = ma / af;
b = 0.0705; %Bore m
b=.1
s = 0.0641; %Stroke m
s=.1
Ao = pi*(b/2)^2; %Total Area of Combustion chamber at TDC, will come from 
Vd = 0.0015/6;
V1 = r/(r-1)*Vd;
Vo = V1/r; %Total volume at top dead center, m^3
Vr = V1; %Volume at the time the intake valve closes
P1 = 101325; %Initial Pressure, Pa
P1 = 100000
Pr = P1; %Pressure when intatke valve closes; assumed that the intake valve closes at bottom dead center. 
T1 = 300; %K
Tr = T1; %K
Tw = 525; %K
Tbar = Tw / T1;
Tbar = 1.2
Qlhv = 44000000; %J/kg
Qin = mf*Qlhv/(P1*V1); %Dimensionless total Heat release
Qin = 27.22
N = 5000; %rpm
omega = N*2*pi/60; %rps
omega = 314.1
c = 0; %mass loss coefficient
Ubar = 2*s*N/60; %Mean piston speed (m/s)
Ubar = s*omega/pi
B = 4*r/(b*(Ao/Vo)-4); %dimensionless volume
B = 2.22

%Engine Geometry & parameters
g(1) = r;
g(2) = Vd;
g(3) = Vr;
g(4) = b;
g(5) = Ubar;
g(6) = B;
g(7) = P1;
g(8) = Pr;
g(9) = Qin;
g(10) = T1;
g(11) = Tr;
g(12) = thetas;
g(13) = thetad;
g(14) = omega;
g(15) = Tbar;

dxb_dtheta = zeros(NN,1);
dV_dtheta = zeros(NN,1);

prop.theta = thetai:step:thetaf;
prop.theta = prop.theta.';
prop.vol = zeros(NN,1);
prop.area = zeros(NN,1);
prop.press = zeros(NN,1);
prop.temp = zeros(NN,1);
prop.rho = zeros(NN,1);
prop.work = zeros(NN,1);
prop.heatloss = zeros(NN,1);
prop.htcoeff = zeros(NN,1);
prop.heatflux = zeros(NN,1);
gam = zeros(NN,1);

gam(1)=calc_gamma(T1);
prop.vol(1) = nondimV(0,r);
prop.area(1) = (Ao-4*Vo/b) + 4* prop.vol(1)/b;
prop.press(1) = 1; %nondimensional initial pressure, P1/P1
prop.temp(1) = T1; %initial temperature, K

for i = 2:NN
    [prop.vol(i),dV_dtheta(i)] = nondimV(prop.theta(i), r);
    prop.rho(i) = ma / (V1*prop.vol(i));
    prop.area(i) = (Ao-4*Vo/b) + 4* prop.vol(i)*V1/b;
end

for i = 2:NN
    gamma = calc_gamma(prop.temp(i-1));
    gamma = 1.4;
    gam(i) = gamma;
    [prop.press(i),prop.heatflux(i-1),prop.htcoeff(i-1)] = rk4(step,prop.theta(i-1),prop.press(i-1),prop.temp(i-1),g,gamma); % Assumes that the change in theta is equal to whatever the unit of theta is. This program is written in degrees, but, if this program were written in radians I think it would work out the same. 
    prop.temp(i) = prop.press(i)*P1/(prop.rho(i)*R);
end

prop.heatflux = prop.heatflux.*T1./(prop.area .*10^6); %MW/m^2

prop.press = prop.press .* P1; %Pa
prop.vol = prop.vol .* V1; %m^3

mprops.press = mprops.press ./ 10^6; %MPa
prop.press = prop.press / 10^6; %MPa
bar = prop.press .* 10; %bar

figure
hold on
plot(prop.theta,prop.press)
plot(mprops.theta,mprops.press)
ylabel("Pressure, MPa")
xlabel("Angle, rad")
title("Pressure vs Crank Angle")
legend("Actual Pressure", "Motored Pressure")
hold off 

figure
plot(prop.theta,bar)
ylabel("Pressure, bar")
xlabel("Angle, rad")
title("Pressure vs Crank Angle")

figure
plot(prop.theta,prop.temp)
ylabel("Temperature, K")
xlabel("Angle, rad")
title("Temperature vs Crank Angle")

figure
plot(prop.theta,gam)
ylabel("Gamma")
xlabel("Angle, rad")
title("Gamma vs Crank Angle")

figure
plot(prop.theta,prop.vol)
ylabel("Volume, m^3")
xlabel("Angle, rad")
title("Volume vs Crank Angle")

figure
plot(prop.theta,prop.heatflux)
ylabel("q'', MW/m^2")
xlabel("Angle, rad")
title("Heat loss vs Crank Angle")

figure
plot(prop.theta,prop.htcoeff)
ylabel("h")
xlabel("Angle, rad")
title("Heat transfer coefficient vs Crank Angle")

%% 
%Approximates motored pressure at the points immediately to the right and
%left as being approximately equivalent to point i
function [P2,dQw_dtheta,ht] = rk4(h,theta,P,T,g,gamma)
    F = [0,0,0,0];
    [F(1), dQw_dtheta, ht] = deriv(theta,P,T,g,gamma);
    F(1) = F(1) * h;
    [F(2), ~,~] = deriv(theta+h/2,P+F(1)/2,T,g,gamma);
    F(2) = F(2)*h;
    [F(3), ~,~] = deriv(theta+h/2,P+F(2)/2,T,g,gamma);
    F(3) = F(3)*h;
    [F(4),~,~] = deriv(theta+h,P+F(3),T,g,gamma);
    F(4) = F(4)*h;
    P2 = P + (F(1)+2*F(2)+2*F(3)+F(4))/6;
end

function [dP_dtheta, dQw_dtheta, h] = deriv(theta,P,T,g,gamma)
    r=g(1); 
    Vd=g(2); 
    Vr=g(3);
    b=g(4); 
    Ubar=g(5); 
    B=g(6); 
    P1=g(7); 
    Pr=g(8); 
    Qin=g(9); 
    T1=g(10);
    Tr=g(11);
    thetas=g(12);
    thetad=g(13);
    w=g(14);
    Tbar=g(15);

    [V,dV_dtheta] = nondimV(theta,r);
    U = 2.28 * Ubar + 0.00324*T1*(r-1)*(P-V^(-gamma))/r;
    if theta < thetas 
        dxb_dtheta = 0;
    else
        [~, dxb_dtheta] = burnRate(theta, thetas, thetad);
    end
    
    h = heatTransferCoeff(P*P1/1000,U,b,T);
    dQw_dtheta = HeatLoss(P,V,B,h,b,T1,P1,Tbar,w);
    
    dP_dtheta = PressureChange(gamma,P,V,Qin,dxb_dtheta,dQw_dtheta,dV_dtheta);
end

function rate = PressureChange(gamma, P, V, Qin, dx_dtheta, dQw_dtheta, dV_dtheta)
rate = (gamma-1)/V*(Qin*dx_dtheta-dQw_dtheta)-gamma*P/V*dV_dtheta;
end

function dQw_dtheta = HeatLoss(P,V,B,h,b,T1,P1,Tbar,w)
    h_nd = 4*h*T1/(P1*w*B*b); %nondimensional heat transfer coefficient, nondimensionalized by values at state 1. 
    dQw_dtheta = h_nd*(1+B*V)*(P*V - Tbar); %heat loss through the cylinder wall
end

function h = heatTransferCoeff(P,U,b,T)
h = 3.26*(P^0.8)*(U^0.8)*(b^(-.2))*(T^(-.55));
end

function [V,rate] = nondimV(theta,r)
    V = (1/r) + (r-1)/(2*r).*(1+cos(theta));
    rate = -(r-1)/(2*r)*sin(theta);
end

function [xb, rate] = burnRate(theta, thetas, thetad)
e = 2.712821282;
a = 5;
n = 3;
xb = 1-e^(-a*((theta - thetas)/thetad)^n);
rate = n*a*(1-xb)*((theta - thetas)/thetad)^(n-1);
end