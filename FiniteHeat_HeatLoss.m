%NOTE: Matlab trig functions take angles in radians as arguments, but this
%program uses angles in degrees. This needs to be converted before this
%program is implemented. 

%Real Engine Cycle with Heat Transfer from the Cylinder
step = 1; %crankangle interval for calculations
NN = 360/step;  %number of datapoints
theta = -180; %initial crankangle
thetae = theta + step; %final crankangle througout steps

thetas = -20; %start of heat release (deg)
thetad = 40; %duration of heat release (deg)
num = thetad/step;
xb = zeros(num,1);
r = 10; %compression ratio
ma = ; %Get from intake/outtake program
af = 14.7; %Air to fuel ratio
mf = ma / af;
Ao = pi*(b/2)^2; %Total Area of Combustion chamber at TDC, will come from 
Vo = 0.00015; %Total volume at top dead center, m^3
V1 = 0.00165; %Total volume at bottom dead center
P1 = 101325; %Initial Pressure, Pa
T1 = 300; %K
Qin = mf*Qlhv/(P1*V1); %Dimensionless total Heat release
a = 5;
n = 3;
N = 5000; %rpm
omega = N/60; %rps
c = 0; %mass loss coefficient
s = 0.0705; %Stroke, m
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
prop.mass = zeros(NN,1);
prop.htcoeff = zeros(NN,1);
prop.heatflux = zeros(NN,1);

for i = 2:NN
    gamma = calc_gamma(prop.temp(i));
    prop.theta(i) = theta;
    prop.vol(i) = nondimV(prop.theta(i), r);
    prop.rho(i) = ma / prop.vol(i);
    prop.area(i) = (Ao-4*Vo/b) + 4* prop.vol(i)/b;
    dV_dtheta(i) = (r-1)/(2*r)*sin(theta);

    if theta < thetas || theta > (thetas + thetad) %VALVES ARE OPEN check this
       U = 6.18*Ubar;
    else %VALVES ARE CLOSED Check this
       U = 2.28 * Ubar + 0.00324*Tr*(Vd/Vr)*(prop.press(i)-___MOTORED PRESSURE AT SAME CRANK ANGLE____) / Pr;
    end
    
    h = 3.26 * (prop.press(i-1) * U)^0.8 * b^(-.2) * prop.temp(i-1)^(-0.55);
    
    if theta < thetas 
        dxb_dtheta(i) = burnRate(prop.theta(i), thetas, thetad,a,n,xb(i));
    elseif theta > (thetad + thetas)
        dxb_dtheta(i) = 1;
    else
        dxb_dtheta(i) = 0;
    end
    
    prop.heatflux(i) = HeatLoss(P,prop.vol(i),B,h,b,T1,P1,Tbar);
    
    dP_dtheta = PressureChange(gamma, prop.press(i-1), prop.vol(i-1), Qin, dx_dtheta(i), dQw_dtheta(i), dV_dtheta(i));
    prop.press(i) = prop.press(i-1) + dP_dtheta * step; % Assumes that the change in theta is equal to whatever the unit of theta is. This program is written in degrees, but, if this program were written in radians I think it would work out the same. 
    prop.temp(i) = prop.press(i)*prop.vol(i)/(prop.rho(i)*R);
end

function rate = PressureChange(gamma, P, V, Qin, dx_dtheta, dQw_dtheta, dV_dtheta)
rate = (gamma - 1)/V * (Qin*dx_dtheta - dQw_dtheta) - gamma * P / V * dV_dtheta;

end

function V = nondimV(theta,r)
    V = (1/r) + (r-1)/(2*r)*(1-cos(theta));
end

function [xb, rate] = burnRate(theta, thetas, thetad,a,n)
e = 2.712821282;
xb = 1-e^(-a*((theta - thetas)/thetad)^n);
rate = n*a*(1-xb)*((theta - thetas)/thetad)^(n-1);
end

function dQw_dtheta = HeatLoss(P,V,B,h,b,T1,P1,Tbar)
    h = 4*h*T1/(P1*w*B*b); %nondimensional heat transfer coefficient, nondimensionalized by values at state 1. 
    dQw_dtheta = h*(1+B*V)*(P*V - Tbar); %heat loss through the cylinder wall
end