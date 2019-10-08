%%
% This program  computes various trade studies, comparing various metrics 
% as a function of displacement and compression ratio
clear all
close all
set(0, 'DefaultAxesFontWeight', 'normal', ...
      'DefaultAxesFontSize', 18, ...
      'DefaultAxesFontAngle', 'normal', ... 
      'DefaultAxesFontWeight', 'normal', ... 
      'DefaultAxesTitleFontWeight', 'bold', ...
      'DefaultAxesTitleFontSizeMultiplier', 1.2) ;
set(groot,'defaultLineLineWidth',3)
i = 1 
figure
for rc = 9:.1:10
    
    % Initial Conditions
    p(1,i) = 101325;%pa
    %p(1,i) = p(1)*1000; %pa
    T(1,i) = 300; %k
    %(1,i) = 300; %k
    gamma = 1.4;
    R = 287; %j/kg, gas constant
    %rc = 9; % compression ration
    Qlhv = 44e6; %j/kg
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
    for D = .0015 : .000050 : .0018
        % Work
        %m^3
        C = 6; % initial guess
        Cv = Cp-R;
        Ws(i, j) = Cv*((T(3, i)-T(2, i)) - (T(4, i)-T(1, i)));
        Ma(i, j) = (D/C)*(rc/(rc-1)) * rho(1, i);
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
        j = j+ 1;
    end
    otto_eff(i) = 1-(1/rc^(gamma-1))
    o_eff(i) = 1-(T(4,i)-T(1,i))./(T(3,i)-T(2,i))
    i = i + 1; 
end

%%

% Thermal efficiency trade study: Compression ratio
figure
plot(9:.1:10, otto_eff)
ylabel("Thermal efficiency")
hold on
yyaxis right
plot(9:.1:10, SFC_Converted(:, 1));
   
%st = strcat('Compression Ratio =  ' , num2str(9 + .1*(i-1)))
%text(.0016,SFC_Converted(i,1)+.25, st);
%plot(SFC_Converted(i,:)
%st = strcat('Compression Ratio =  ' , num2str(9.0))
%text(.0016,SFC_Converted(1,1)+.25, st);
ylabel('Specific Fuel Consumption (g/Kw-hr)')
xlabel('RC')
title('Specific Fuel Consumption and Thermal Efficiency vs Compression Ratio')

rho_1 = p(1,1)/(R*T(1,1))
mass_flow_air = SFC(1,1)*Power(1,1)*f

vol_eff = 2*mass_flow_air/(rho_1*.0015*(N/60))

%%

% Trade Study of compression ratio and displacement vs horsepower
figure
for i = 1 : 7
    plot(P_hp(:,i),9:.1:10);
    hold on
    st = strcat('Displacement =  ' , num2str(.00150 + .00005*(i-1)), ' cc');
    ht = text(P_hp(1, i ),9.3, st);
    set(ht, 'Rotation', 85);
end
xlabel('Horse Power')
ylabel('Compression ratio')
title('Horse Power vs Compression Ratio and Displacement')
%%
C = 6
air_flow = Ma * N * C/120

figure 
hold on

for i = 1 : 7
    plot(air_flow(:,i),9:.1:10);
    st = strcat('Displacement =  ' , num2str(.00150 + .00005*(i-1)), ' cc');
    ht = text(air_flow(1, i ),9.3, st);
    set(ht, 'Rotation', 95);
end

xlabel('Air flow (kg/s)')
ylabel('Compression ratio')
title('Airflow vs Compression Ratio and Displacement')
xlim([.080, .1])