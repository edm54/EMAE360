%%
%clear all
% Requires a run of valve_lift before running for curves
close all
clear 

set(0, 'DefaultAxesFontWeight', 'normal', ...
      'DefaultAxesFontSize', 18, ...
      'DefaultAxesFontAngle', 'normal', ... 
      'DefaultAxesFontWeight', 'normal', ... 
      'DefaultAxesTitleFontWeight', 'bold', ...
      'DefaultAxesTitleFontSizeMultiplier', 1.2) ;
set(groot,'defaultLineLineWidth',3)

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


figure
equivolence_ratio = 1.1

    % Equivolence Ratio * Ideal Air Fuel 
    fuel_air = 14.7/equivolence_ratio

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

    combustion_eff = combustionEff(equivolence_ratio)
    T(3) = 2868.97; %K, from CEA
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
    stroke = .0641;
    bore = 1.1* stroke; 
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

    % Equivolence Ratio * Ideal Air Fuel 
    fuel_air = 14.7/equivolence_ratio

    C = 6; % number of cylinders
    D = .0015; %Total Displacement, m^3

    % State 1, 
    rho(1) = p(1)/(R*T(1)); %from ideal gas law

    % combustive efficiency, from Heywood p 82
    if f <= 14.7
        nc = 0.95;
    end
    
    N = 800:100:15000;
    Veff = volumetric_efficiency(N);
    i = 1
    %for N = 800:100:10000

    for i = 1:length(N)
        Ws = Cv*((T(3)-T(2)) - (T(4)-T(1))); %Specific work per cylinder, J/kg
        Ma(i) = Veff(i) * (D/C)*(rc/(rc-1)) * rho(1); %Mass of Air in each cylinder, kg
        %Ma(k,i) = (D/C)*(rc/(rc-1)) * rho(1);
        D_cc = 1500;

        %Mf = Ma/f;
        Mf(i) = Ma(i)/fuel_air;
        Qin = combustion_eff * Mf(i) * Qlhv; 
        Vd = ((pi*bore^2*stroke)/4); %meter cubes
        V1  = Vd/(1-(1/rc));
        Q(i) =  Qin /(p(1)*V1);
        net_work(i) = 6 * FiniteHeatRelease(Q(i), N(i), Ma(i), 0);

        %%
        Wt(i) = net_work(i) * (p(1)*V1);
        Power(i) = Wt(i) * N(i)/120;
        P_rate(i) = .8 * Power(i); 
        %efficiency(i) = combustion_eff * V_eff(i) * mechanical_eff;
        efficiency(i) = 1;
        RPM = [2100 15000];
        %mech_eff = [ .9 .75];  
        %mech_eff = [ .8 .75];  
        mech_eff = [ .9 .65];  
        if N(i)<= 2100
            mechanical_eff = .9;
        %elseif N(i)> 10000
           % mechanical_eff = .75
        else
            mechanical_eff = interp1(RPM, mech_eff, N(i), 'linear');
        end
        
        
        efficiency(i) = mechanical_eff;
%         P_specific(i) = efficiency(i) * Wt(k,i) * N(i)/120;
%         P_cylinder(i) = P_specific(i) * Ma(k,i);
        %P_total(i) = P_cylinder(i);
        imep(i) = Power(i)/1000*2*10^3/(D*1000*(N(i)/60));
        bmep(i) =  mechanical_eff*imep(i);
        %p_real(i) = ((P_total(i) - Pf(i))/(P_total(i))) * P_total(i)

        %P_hp(k,i) = P_total(i)/745.699872;
       
        P_hp(i) = efficiency(i) * Power(i)/745.699872;
        % SFC
        SFC(i) = (C * Ma(i)/f)/(Wt(i)); %kg/kj
        SFC_Converted(i) = SFC(i) * 3.6e9; %g/Kw-hr
        Torque(i) = Power(i)/((2*3.1415*N(i)/60));

   
     end


    %hold on
    %plot(N,P_hp(k, :))
   
    %legend('.8', '.9', '1', '1.1', '1.2', '1.3')
%%
figure 

disp(N)
plot(N, P_hp)

title('Lift Power vs RPM')
ylabel('Power (hp)')
xlabel('RPM')
xlim([0 15000])

%%
figure 
hold on

plot(N, Torque(1,:))
 
title('Torque vs RPM')
xlabel('Torque')
ylabel('RPM')


%%
figure
plot(N, Veff)
title('Volumetric Efficiency vs RPM')
xlabel('RPM')
ylabel('Volumetric Efficiency')




