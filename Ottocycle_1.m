%%
% This program calculates the airflow and fuel flow rates as a function of
% RPM based on an ideal ottocycle

clear all
i = 1 
figure
rc = 10
    
% Initial Conditions
p(1, i) = 101325;%pa
%p(1,i) = p(1)*1000; %pa
T(1,i) = 294; %k
%(1,i) = 300; %k
gamma = 1.4;
R = 287; %j/kg, gas constant
%rc = 9; % compression ration
Qlhv = 42.7e6; %j/kg
Cp = R*gamma/(gamma-1); 
f = 14.7; % air to fuel ratio

% State 1
rho(1, i) = p(1, i)/(R*T(1, i)); 

% State 2
rho(2, i) = rho(1, i) * rc;
p(2, i) = p(1, i)* (rho(2, i)/rho(1, i))^gamma;
T(2, i) = p(2, i)/(rho(2, i)* R);

% State 3
T(3, i) = T(2, i) + .8*Qlhv/(Cp* f); %.8 is thermal efficeny of combustionn
rho(3, i) = rho(2, i);
p(3, i) = rho(3, i) * R* T(3,i);

% State 4
rho(4, i) = rho(1, i);
p(4, i) = p(3, i)* (rho(4, i)/rho(3, i))^ (gamma); %this is an estimate, could use refprop to find actual
T(4, i) = p(4, i)/(rho(4, i) * R);
j = 1;
D = .0015
% Work
%m^3
C = 6; % initial guess
Cv = Cp-R;
Ws(i, j) = Cv*((T(3, i)-T(2, i)) - (T(4, i)-T(1, i)));
Ma(i, j) = (D/C)*(rc/(rc-1)) * rho(1, i)
Wc(i, j) = Ws(i, j)*Ma(i, j);
Wt(i, j) = Wc(i, j) * C; %total work

%Power
N = 5000; % RPM, Max
Power(i, j) = Wt(i, j) * N/120;
P_rate(i, j) = .8 * Power(i, j); 
P_specific(i, j) = Ws(i, j) * N/120;
P_cylinder(i, j) = P_specific(i, j) * Ma(i, j);
P_total(i, j) = P_cylinder(i, j) * C;
P_hp(i, j) = P_total(i, j)/745.699872;

% SFC
SFC(i, j) = (C * Ma(i, j)/f)/(Wt(i, j)); %kg/kj
SFC_Converted(i, j) = SFC(i,j) * 3.6e9; %g/Kw-hr

otto_eff(i) = 1-(1/rc^(gamma-1))
o_eff(i) = 1-(T(4,i)-T(1,i))./(T(3,i)-T(2,i))

%%
figure
C = 6
air_flow = Ma * N * C/120

i = 1
for N = 800: 100 : 9000
    air_flow_rpm(i) = Ma * N * C/120
    fuel_rate(i) = air_flow_rpm(i)/f
    i = i + 1
end

plot(800: 100 : 9000, air_flow_rpm)
hold on
plot(800: 100 : 9000, fuel_rate)
title('Air flow & Fuel flow vs RPM')
xlabel('RPM')
ylabel('Mass flow rate, air and fuel (kg/s)')
legend('Air mass flow rate', 'Fuel mass flow rate')

