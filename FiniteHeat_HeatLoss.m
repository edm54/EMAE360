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
Vo = 0.00015; %Total volume at top dead center, m^3
V1 = 0.00165; %Total volume at bottom dead center
Qin = ; %Dimensionless total Heat release
a = 5;
n = 3;
N = 5000; %rpm
omega = N/60; %rps
c = 0; %mass loss coefficient
s = 0.0705; %Stroke, m
b = 0.0641; %Bore, m
T_bdc = 300; %Temperature at bdc (K)
P_bdc = 100 ;%Pressure at bdc (K)
Up = s*omega/pi; %Mean piston speed (m/s)
beta = 4*Vo/(b*(Ao-4*Vo/b)); %dimensionless volume

properties.theta = zeros(NN,1);
properties.vol = zeros(NN,1);
properties.press = zeros(NN,1);
properties.work = zeros(NN,1);
properties.heatloss = zeros(NN,1);
properties.mass = zeros(NN,1);
properties.htcoeff = zeros(NN,1);
properties.heatflux = zeros(NN,1);

for i = 1:NN
    properties.theta(1) = theta;
    properties.vol(1) = nondimV(theta);
    dxb_dtheta = burnRate(properties.theta(i), thetas, thetad,a,n,xb(i));
    
    dP_dtheta = pressureChange();
end

function rate = PressureChange()
gamma = calcgamma(T);
rate = (gamma - 1)/V * (Qin*dx_dtheta - dQw_dtheta) - gamma * P / V * dV_dtheta;

end

function V = nondimV(theta)

end

function rate = burnRate(theta, thetas, thetad,a,n,xb)
rate = n*a*(1-xb)*((theta - thetas)/thetad)^(n-1);
end