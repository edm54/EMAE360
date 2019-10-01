clear all
close all
i = 1; 

rc = 10;
    
% Initial Conditions
p(1) = 101325;%pa
T(1) = 294; %k
gamma = 1.4;
R = 287; %j/kg*K, gas constant
Cp = R*gamma/(gamma-1); %J/kg*K
Cv = Cp-R;

%Engine Parameters
Qlhv = 42.7e6; %j/kg
f = 14.7; % air to fuel ratio
C = 6; % number of cylinders

% State 1, 
rho(1) = p(1)/(R*T(1)); %from ideal gas law
v(1) = 1/rho(1);

% State 2, path 2-3: isentropic compression
rho(2) = rho(1) * rc;
p(2) = p(1)* (rho(2)/rho(1))^gamma;
T(2) = p(2)/(rho(2)* R);
u_su = 390; %kJ/kg, from figure 4.3, T = 740K
v(2) = 1/rho(2);

%Combustion Calculations
equiv = 1;
xb = 0.02;
v(3) = v(2);
u_form = -118.2 - 29568*xb;
h_form = -129.7-2951*xb;
ub = u_su + u_form;
Tb = 2700;
pb = 10000; %kPa