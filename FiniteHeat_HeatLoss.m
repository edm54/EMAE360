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
mstep = 0.0005;
thetaf = 2*pi;
mprops.theta = thetai:mstep:thetaf;
mprops.theta = mprops.theta.';
mprops.vol = zeros(length(mprops.theta),1);
mprops.press = zeros(length(mprops.theta),1);
mprops.temp = zeros(length(mprops.theta),1);
mprops.rho = zeros(length(mprops.theta),1);
gamma = zeros(length(mprops.theta),1);

mprops.vol(1) = nondimV(thetai, r);
mprops.press(1) = 1;
mprops.temp(1) = T1;
k = P1*(V1)^(calc_gamma(T1));
mprops.rho(1) = 1.2;
ma = mprops.vol(1)*V1/mprops.rho(1);

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

xb = zeros(length(num),1);

r = 10; %compression ratio
ma = 0.000336; %kg, Get from intake/outtake program
af = 14.7; %Air to fuel ratio
mf = ma / af;
s = 0.0705; %Stroke m
b = 0.0641; %Bore m
Ao = pi*(b/2)^2; %Total Area of Combustion chamber at TDC, will come from 
Vo = 0.00015; %Total volume at top dead center, m^3
V1 = 0.00165; %Total volume at bottom dead center
Vr = V1; %Volume at the time the intake valve closes
Vd = 0.0015; %Displacement Volume
P1 = 101325; %Initial Pressure, Pa
Pr = 101325; %Pressure when intatke valve closes; assumed that the intake valve closes at bottom dead center. 
T1 = 300; %K
Tr = T1; %K
Tw = 550; %K
Tbar = Tw / T1;
Qlhv = 44000000; %J/kg
Qin = mf*Qlhv/(P1*V1); %Dimensionless total Heat release
N = 5000; %rpm
omega = N/60; %rps
c = 0; %mass loss coefficient
Ubar = s*omega/pi; %Mean piston speed (m/s)
B = 4*Vo/(b*(Ao-4*Vo/b)); %dimensionless volume

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

prop.vol(1) = nondimV(thetas,r);
prop.area(1) = (Ao-4*Vo/b) + 4* prop.vol(1)/b;
prop.press(1) = 1; %nondimensional initial pressure, P1/P1
prop.temp(1) = T1; %nondimensional initial temperature, T1/T1
prop.heatflux = zeros(0);

for i = 2:NN
    [prop.vol(i),dV_dtheta(i)] = nondimV(prop.theta(i), r);
    prop.rho(i) = ma / (V1*prop.vol(i));
    prop.area(i) = (Ao-4*Vo/b) + 4* prop.vol(i)/b;
end

for i = 2:NN
    gamma = calc_gamma(prop.temp(i-1)*T1);
    prop.press(i) = rk4(step,prop.theta(i-1),prop.press(i-1),mprops.press(2*i-1),mprops.press(2*i),mprops.press(2*i+1),prop.temp(i-1),g,gamma); % Assumes that the change in theta is equal to whatever the unit of theta is. This program is written in degrees, but, if this program were written in radians I think it would work out the same. 
    prop.temp(i) = prop.press(i)*P1*V1*prop.vol(i)/(prop.rho(i)*R);
end

%% 
function P2 = rk4(h,theta,P,Pm,Pm1,Pm2,T,g,gamma)
    F1 = h*deriv(theta,P,Pm,T,g,gamma);
    F2 = h*deriv(theta+h/2,P+F1/2,Pm1,T,g,gamma);
    F3 = h*deriv(theta+h/2,P+F2/2,Pm1,T,g,gamma);
    F4 = h*deriv(theta+h,P+F3,Pm2,T,g,gamma);
    P2 = P + (F1+2*F2+2*F3+F4)/6;
end

function dP_dtheta = deriv(theta,P,Pm,T,g,gamma)
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

    U = 2.28 * Ubar + 0.00324*Tr*(Vd/Vr)*(P*P1-Pm) / Pr;
    if theta < thetas 
        dxb_dtheta = 0;
    elseif theta > (thetad + thetas)
        dxb_dtheta = 1;
    else
        [~, dxb_dtheta] = burnRate(theta, thetas, thetad);
    end
    
    [V,dV_dtheta] = nondimV(theta,r);
    
    h = heatTransferCoeff(P*P1,U,b,T);
    dQw_dtheta = HeatLoss(P,V,B,h,b,T1,P1,Tbar,w);
    
    dP_dtheta = PressureChange(gamma,P,V,Qin,dxb_dtheta,dQw_dtheta,dV_dtheta);
end

function rate = PressureChange(gamma, P, V, Qin, dx_dtheta, dQw_dtheta, dV_dtheta)
rate = (gamma - 1)/V * (Qin*dx_dtheta - dQw_dtheta) - gamma * P / V * dV_dtheta;

end

function dQw_dtheta = HeatLoss(P,V,B,h,b,T1,P1,Tbar,w)
    h_nd = 4*h*T1/(P1*w*B*b); %nondimensional heat transfer coefficient, nondimensionalized by values at state 1. 
    dQw_dtheta = h_nd*(1+B*V)*(P*V - Tbar); %heat loss through the cylinder wall
end

function h = heatTransferCoeff(P,U,b,T)
h = 3.26*P^0.8*U^0.8*b^(-.2)*T^(-.55);
end

function [V,rate] = nondimV(theta,r)
    V = (1/r) + (r-1)/(2*r).*(1+cos(theta));
    rate = (r-1)/(2*r)*sin(theta);
end

function [xb, rate] = burnRate(theta, thetas, thetad)
e = 2.712821282;
a = 5;
n = 3;
xb = 1-e^(-a*((theta - thetas)/thetad)^n);
rate = n*a*(1-xb)*((theta - thetas)/thetad)^(n-1);
end