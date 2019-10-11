%%
clear all
close all
rc = 10;
    
% Initial Conditions
p(1) = 101325;%pa
T(1) = 294; %k
gamma = calc_gamma(T(1));
R = 287; %j/kg*K, gas constant
Cp = R*gamma/(gamma-1); %J/kg*K
Cv = Cp-R;

%Engine Parameters
Qlhv = 44e6; %j/kg from Heywood App D
f = 14.7; % air to fuel mass ratio
C = 6; % number of cylinders
D = .0015; %Total Displacement, m^3

% State 1, 
rho(1) = p(1)/(R*T(1)); %from ideal gas law

% State 2, path 2-3: isentropic compression
rho(2) = rho(1) * rc;

interval = (rho(2)-rho(1))/500;
rho_comp = rho(1):interval:rho(2);
pcomp = zeros(1, length(rho_comp));
Tcomp = zeros(1, length(rho_comp));

pcomp(1) = p(1);
Tcomp(1) = T(1);
for i = 2:length(rho_comp)
    gamma = calc_gamma(Tcomp(i-1));
    pcomp(i) = pcomp(i-1)*(rho_comp(i)/rho_comp(i-1))^(gamma);
    Tcomp(i) = pcomp(i) / (rho_comp(i)*R);
end

p(2) = pcomp(length(pcomp));
T(2) = Tcomp(length(Tcomp));

% State 3, path 2-3: constant volume heat addition

% combustive efficiency, from Heywood p 82
if f <= 14.7
    nc = 0.95;
end


T(3) = 2869.; %K, from CEA
rho(3) = rho(2);
p(3) = rho(3) * R* T(3);


% State 4
rho(4) = rho(1);

interval = (rho(4)-rho(3))/500;
rho_comp = rho(3):interval:rho(4);
pcomp = zeros(1, length(rho_comp));
Tcomp = zeros(1, length(rho_comp));

pcomp(1) = p(3);
Tcomp(1) = T(3);
for i = 2:length(rho_comp)
    gamma = calc_gamma(Tcomp(i-1));
    pcomp(i) = pcomp(i-1)*(rho_comp(i)/rho_comp(i-1))^(gamma);
    Tcomp(i) = pcomp(i) / (rho_comp(i)*R);
end

p(4) = pcomp(length(pcomp)); %this is an estimate, could use refprop to find actual
T(4) = Tcomp(length(Tcomp));
gamma = calc_gamma(T(4));

j = 1;
% Work

Ws = Cv*((T(3)-T(2)) - (T(4)-T(1))); %Specific work per cylinder, J/kg
Ma = (D/C)*(rc/(rc-1)) * rho(1); %Mass of Air in each cylinder, kg
Wc = Ws*Ma; %Work per Cylinder, J
Wt = Wc * C; %total work, J
nv = Ma*C/(rho(1)*D);
nt = 1-(T(4)-T(1))/(T(3)-T(2));

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
i = 1;
%% Mechanical eff

stroke = stroke / 100;
RPM = [2100 9000];
mech_eff = [ .9 .75];
%maximum rpm
Nmax = 7000;

stroke = stroke / 100;
for N = 1500:25:9400
    c = 6; % number of cylinder
%{
    Ubar = 2*(stroke) * N/60;

    A = pi * bore * stroke/100;

    m = .17; % mm
    m = .07/1000; %m
    stroke = stroke / 100;


    %u = (8.34* 10^-5) * e^(1474/(temp - 368)
    u = 1.44E-4;
    %https://wiki.anton-paar.com/en/engine-oil/
    u = .0119;

    t = u * (Ubar / m);
    f = A * t;

    Pf(i) = 1.5 * c * Ubar * f;
    %}
    RPM = [2100 9400];
    mech_eff = [ .9 .75];   
    if N<= 2100
        mechanical_eff = .9;
    else
        mechanical_eff = interp1(RPM, mech_eff, N, 'linear');
    end
    Power(i) = Wt * N/120;
    P_rate(i) = .8 * Power(i); 
    P_specific(i) = mechanical_eff* Ws * N/120;
    P_cylinder(i) = P_specific(i) * Ma;
    P_total(i) = P_cylinder(i) * C;
    imep(i) = P_total(i)/1000*2*10^3/(D*1000*(N/60));
    bmep(i) = mechanical_eff*imep(i);
    %p_real(i) = ((P_total(i) - Pf(i))/(P_total(i))) * P_total(i)
    
    P_hp(i) = P_total(i)/745.699872;
    % SFC
    SFC(i) = (C * Ma/f)/(Wt); %kg/kj
    SFC_Converted(i) = SFC(i) * 3.6e9; %g/Kw-hr
    Torque(i) = P_total(i)/((2*3.1415*N/60));
    
    i = i + 1;
end

fprintf('                  Mountain CWRU: Baja Blast           \n');
fprintf(' Outputs at 5000 RPM \n')
fprintf(' Net Work                 %7.2f             \n', Wt);
fprintf(' Net Power                %7.2f \n', P_hp(141))
fprintf(' Torque         %5.3f               \n', Torque(141));
fprintf(' Thermal Efficiency         %5.3f               \n', nt);
fprintf(' Bmep                      %5.2f                \n', bmep(141));
fprintf(' Outputs at 9400 RPM \n')
fprintf(' Power             %5.2f \n', P_hp(length(P_hp)))

pl = input('Plot? Type y or n \n ->', 's');

if strcmp(pl, 'y')
    figure 
    plot(1500:25:9400, P_total);
    xlabel('RPM')
    ylabel('Total power, Watts')
    title('Power as a function of RPM')

    figure 
    P_spec_hp = P_specific/745.1;
    plot(1500: 25:9400, P_spec_hp);
    xlabel('RPM')
    ylabel('Specific power, hp/m^2')
    title('Specific power as a function of RPM')


    figure 
    plot(1500: 25:9400, Torque);
    xlabel('RPM')
    ylabel('Torque, N*m')
    title('Torque vs RPM')


    %{
    figure 
    plot(1500: 25:7500, Pf);
    xlabel('RPM')
    ylabel('Friction loss')
    title('Friction loss vs. RPM')
    %S(1, i)= refpropm('S','D',rho(1,i), 'P',p(1, i)/1e3, 'air.ppf');
    %S(2, i)= refpropm('S','D',rho(2, i), 'P',p(2, i)/1e3, 'air.ppf');
    %S(3, i)= refpropm('S','D',rho(3, i), 'P',p(3, i)/1e3, 'air.ppf');
    %S(4, i)= refpropm('S','D',rho(4, i), 'P',p(4, i)/1e3, 'air.ppf');
    %}

    figure
    plot(1500:25:9400, P_hp)
    xlabel('RPM')
    ylabel('Power (hp)')
    title('Power vs RPM')

    figure
    plot(1500:25:9400, bmep)
    xlabel('RPM')
    ylabel('bmep (kPa)')
    title('bmep vs RPM')
end
otto_eff = 1-(1/rc^(gamma-1));
o_eff = 1-(T(4)-T(1))./(T(3)-T(2));
 

