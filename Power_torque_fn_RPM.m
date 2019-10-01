%%
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

% State 2, path 2-3: isentropic compression
rho(2) = rho(1) * rc;
p(2) = p(1)* (rho(2)/rho(1))^gamma;
T(2) = p(2)/(rho(2)* R);

% State 3, path 2-3: constant volume heat addition

% combustive efficiency, from Heywood p 82
if f <= 14.7
    nc = 0.95;
end

T(3) = T(2) + nc*Qlhv/(Cp* f); 
rho(3) = rho(2);
p(3) = rho(3) * R* T(3);

% State 4
rho(4) = rho(1);
p(4) = p(3)* (rho(4)/rho(3))^ (gamma); %this is an estimate, could use refprop to find actual
T(4) = p(4)/(rho(4) * R);
j = 1;
D = .0015; %Total Displacement, m^3
% Work

Ws = Cv*((T(3)-T(2)) - (T(4)-T(1))); %Specific work per cylinder, J/kg
Ma = (D/C)*(rc/(rc-1)) * rho(1); %Mass of Air in each cylinder, kg
Wc = Ws*Ma; %Work per Cylinder, J
Wt = Wc * C; %total work, J
nv = Ma*C/(rho(1)*D);

RPM = 5000;
D_cc = 1500;

Mf = Ma/f;
Qin = Mf * Qlhv;
S = 2000; %cM/sec
N = RPM/60; %rev/sec
stroke = S/(2 * N); %cM
displacement = D_cc/6;
bore = sqrt(4 * displacement /(pi*stroke));

stroke = .0641;
bore = 1.1* stroke; 

Vd = ((pi*bore^2*stroke)/4); %meter cubes
V1  = Vd/(1-(1/rc));
Q = Qin /(p(1)*V1);

%Power
N = 5000; % RPM, Max
i = 1
%% Mechanical eff

stroke = stroke / 100;
RPM = [2100 9000];
mech_eff = [ .9 .75];
%maximum rpm
Nmax = 7000;

%mechanical efficiency assuming linear variation between Nmax and N = 2100
%rpm
nm = 0.9 + (0.75-0.9)/(Nmax-2100)*(N-2100);

stroke = stroke / 100;
for N = 1500:25:7500
    c = 6; % number of cylinder
   
    Ubar = 2*(stroke) * N/60;

    A = pi * bore * stroke/100;

    m = .17; % mm
    m = .07/1000; %m
    stroke = stroke / 100
    RPM = [2100 9000]
    mech_eff = [ .9 .75]

    %u = (8.34* 10^-5) * e^(1474/(temp - 368)
    u = 1.44E-4;
    %https://wiki.anton-paar.com/en/engine-oil/
    u = .0119;

    t = u * (Ubar / m);
    f = A * t;

    Pf(i) = 1.5 * c * Ubar * f
    
    if N<= 2100
        mechanical_eff = .9
    else
        mechanical_eff = interp1(RPM, mech_eff, N, 'linear')
    end
    Power(i) = Wt * N/120;
    P_rate(i) = .8 * Power(i); 
    P_specific(i) = mechanical_eff* Ws * N/120;
    P_cylinder(i) = P_specific(i) * Ma;
    P_total(i) = P_cylinder(i) * C;
    %p_real(i) = ((P_total(i) - Pf(i))/(P_total(i))) * P_total(i)
    
    P_hp(i) = P_total(i)/745.699872;
    % SFC
    SFC(i) = (C * Ma/f)/(Wt); %kg/kj
    SFC_Converted(i) = SFC(i) * 3.6e9; %g/Kw-hr
    Torque(i) = 9548.8 * P_total(i)/N;
    
    i = i + 1;
end

figure 
plot(1500: 25:8500, P_total);
xlabel('RPM')
ylabel('Total power, Watts')
title('Power as a function of RPM')

figure 
plot(1500: 25:8500, P_specific);
xlabel('RPM')
ylabel('Specific power, Watts/M^2')
title('Specific power as a function of RPM')


figure 
plot(1500: 25:8500, Torque);
xlabel('RPM')
ylabel('Torque, NM')

figure 
plot(1500: 25:8500, Pf);
xlabel('RPM')
ylabel('Friction loss')
title('Friction loss vs. RPM')
%S(1, i)= refpropm('S','D',rho(1,i), 'P',p(1, i)/1e3, 'air.ppf');
%S(2, i)= refpropm('S','D',rho(2, i), 'P',p(2, i)/1e3, 'air.ppf');
%S(3, i)= refpropm('S','D',rho(3, i), 'P',p(3, i)/1e3, 'air.ppf');
%S(4, i)= refpropm('S','D',rho(4, i), 'P',p(4, i)/1e3, 'air.ppf');


%plot(P_hp(i,:),.0015 : .000050 : .0018 )
otto_eff = 1-(1/rc^(gamma-1));
o_eff = 1-(T(4)-T(1))./(T(3)-T(2));
 
%st = strcat('Displacement =  ' , num2str(.00150 + .00005*(i-1)), ' cc')
%ht = text(SFC_Converted(1, i ),9.3, st);
%set(ht, 'Rotation', 85)  
