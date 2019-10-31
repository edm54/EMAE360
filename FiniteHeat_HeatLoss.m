%Motored Engine
clear
close all
r = 10;
T1 = 300; %K 
P1 = 101325; %Pa
V1 = 0.001650; %m^3
ma = 0.000336; %kg
R = 287; %J/kg K

thetai = 0;
step = 0.0005;
thetaf = 2*pi;
mprops.theta = thetai:step:thetaf;
mprops.theta = mprops.theta.';
mprops.vol = zeros(length(mprops.theta),1);
mprops.press = zeros(length(mprops.theta),1);
mprops.temp = zeros(length(mprops.theta),1);
mprops.rho = zeros(length(mprops.theta),1);
gamma = zeros(length(mprops.theta),1);
P = zeros(length(mprops.theta),1);

mprops.vol(1) = nondimV(thetai, r);
mprops.press(1) = 1;
mprops.temp(1) = 1;
k = P1*(V1)^(calc_gamma(T1));
mprops.rho(1) = 1.2;
ma = mprops.vol(1)*V1/mprops.rho(1);

for i = 2:length(props.theta)
    gamma(i) = calc_gamma(mprops.temp(i-1)*T1);
    mprops.vol(i) = nondimV(mprops.theta(i),r) * V1;
    mprops.rho(i) = ma / (mprops.vol(i));
    mprops.press(i) = k / (mprops.vol(i))^(gamma(i));
    mprops.press(i) = mprops.press(i)*P1;
    mprops.temp(i) = mprops.press(i)/(mprops.rho(i)*R);
end

%Real Engine Cycle with Heat Transfer from the Cylinder


%%Overall heat Transfer coefficient

%% Real Engine Cycle with Heat Transfer from the Cylinde
pi = 3.142;
step = 0.001; %crankangle interval for calculation
theta = 0; %initial crankangle, rad
thetae =2*pi; %final crankangle througout steps
NN = (thetae - theta)/step;

thetas = (8/9)*pi; %start of heat release (deg)
thetad = pi/3; %duration of heat release (deg)
num = thetas:step:(thetas + thetad);

xb = zeros(length(num),1);
r = 10; %compression ratio
ma = 0.000336; %kg, Get from intake/outtake program
af = 14.7; %Air to fuel ratio
mf = ma / af;
Ao = pi*(b/2)^2; %Total Area of Combustion chamber at TDC, will come from 
Vo = 0.00015; %Total volume at top dead center, m^3
V1 = 0.00165; %Total volume at bottom dead center
P1 = 101325; %Initial Pressure, Pa
Pr = 101325; %Pressure when intatke valve closes; assumed that the intake valve closes at bottom dead center. 
T1 = 300; %K
Qlhv = 44000000; %J/kg
Qin = mf*Qlhv/(P1*V1); %Dimensionless total Heat release
N = 5000; %rpm
omega = N/60; %rps
c = 0; %mass loss coefficient
s = 0.0705; %Stroke, m+ 
b = 0.0641; %Bore, m
Ubar = s*omega/pi; %Mean piston speed (m/s)
B = 4*Vo/(b*(Ao-4*Vo/b)); %dimensionless volume

dxb_dtheta = zeros(NN,1);
dV_dtheta = zeros(NN,1);

prop.theta = zeros(NN,1);
prop.vol = zeros(NN,1);
prop.area = zeros(NN,1);
prop.press = zeros(NN,1);
prop.temp = zeros(NN,1);
prop.rho = zeros(NN,1);
prop.work = zeros(NN,1);
prop.heatloss = zeros(NN,1);
prop.htcoeff = zeros(NN,1);
prop.heatflux = zeros(NN,1);

prop.theta(1) = thetas;
prop.vol(1) = nondimV(thetas,r);
prop.area(1) = (Ao-4*Vo/b) + 4* prop.vol(1)/b;
prop.press(1) = 1; %nondimensional initial pressure, P1/P1
prop.temp(1) = 1; %nondimensional initial temperature, T1/T1
prop.rho(1) = ma/prop.vol(1);
prop.heatflux = zeros(0);

for i = 2:NN
    prop.theta(i) = prop.theta(i-1) + step;
    prop.vol(i) = nondimV(prop.theta(i), r);
    dV_dtheta(i) = nonDimdV_dtheta(r, theta);
    prop.rho(i) = ma / (V1*prop.vol(i));
    prop.area(i) = (Ao-4*Vo/b) + 4* prop.vol(i)/b;
end

for i = 2:NN
    gamma = calc_gamma(prop.temp(i-1));

    U = 2.28 * Ubar + 0.00324*Tr*(Vd/Vr)*(prop.press(i)-mprop.press(i)) / Pr;
    prop.htcoeff(i) = heatTransferCoeff(prop.press(i-1),U,b,prop.temp(i-1));
       
    prop.heatflux(i) = HeatLoss(P,prop.vol(i),B,prop.htcoeff(i),b,T1,P1,Tbar);
    
    dP_dtheta = PressureChange(gamma, prop.press(i-1), prop.vol(i-1), Qin, dx_dtheta(i), dQw_dtheta(i), dV_dtheta(i));
    prop.press(i) = prop.press(i-1) + dP_dtheta * step; % Assumes that the change in theta is equal to whatever the unit of theta is. This program is written in degrees, but, if this program were written in radians I think it would work out the same. 
    prop.temp(i) = prop.press(i)*prop.vol(i)/(prop.rho(i)*R);
end

function rate = PressureChange(gamma, P, V, Qin, dx_dtheta, dQw_dtheta, dV_dtheta)
rate = (gamma - 1)/V * (Qin*dx_dtheta - dQw_dtheta) - gamma * P / V * dV_dtheta;

end

function V = nondimV(theta,r)
    V = (1/r) + (r-1)/(2*r).*(1+cos(theta));
end

function dV_dT = nonDimdV_dtheta(r, theta)
    dV_dtheta = (r-1)/(2*r)*sin(theta);
end

function [xb, rate] = burnRate(theta, thetas, thetad)
e = 2.712821282;
a = 5;
n = 3;
xb = 1-e^(-a*((theta - thetas)/thetad)^n);
rate = n*a*(1-xb)*((theta - thetas)/thetad)^(n-1);
end

function dQw_dtheta = HeatLoss(P,V,B,h,b,T1,P1,Tbar)
    h_nd = 4*h*T1/(P1*w*B*b); %nondimensional heat transfer coefficient, nondimensionalized by values at state 1. 
    dQw_dtheta = h_nd*(1+B*V)*(P*V - Tbar); %heat loss through the cylinder wall
end

function h = heatTransferCoeff(P,U,b,T)
h = 3.26*P^0.8*U^0.8*b^(-.2)*T^(-.55);
end

function P = rk4(h,
    [xb, dxb_dtheta] = burnRate(theta, thetas, thetad);
    
    if theta < thetas 
        dxb_dtheta(i) = 0;
    elseif theta > (thetad + thetas)
        dxb_dtheta(i) = 1;
    else
        dxb_dtheta(i) = burnRate(prop.theta(i), thetas, thetad,a,n,xb(i));
    end
end