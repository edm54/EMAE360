%%
%clear all
% Requires a run of valve_lift before running for curves
close all
%clear 
P_hp1 = P_hp

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
    %equivolence_ratio = 1.1
    equivolence_ratio = 1.1
    % Equivolence Ratio * Ideal Air Fuel 
    fuel_air = 14.98/equivolence_ratio

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
    
    C = 6; % number of cylinders
    D = .0015; %Total Displacement, m^3

    % State 1, 
    rho(1) = p(1)/(R*T(1)); %from ideal gas law
    
    N = 800 : 100 :9400
    %N = 5000
    Ma = zeros(1, length(N));
    Mf = zeros(1, length(N));
    Q = zeros(1, length(N));
    net_work = zeros(1, length(N));
    Wt = zeros(1, length(N));
    Power = zeros(1,length(N));
    efficiency = zeros(1,length(N));
    imep = zeros(1,length(N));
    bmep = zeros(1,length(N));
    P_hp = zeros(1,length(N));
    SFC = zeros(1,length(N));
    SFC_Converted = zeros(1,length(N));
    Torque = zeros(1,length(N));
    [Veff, pump_work] = volumetric_efficiency(N);
    i = 1;
    %for N = 800:100:10000
    velocity_arr = [30 60 150];
    velocity_arr = [30  150];
    velocity_kmhr = velocity_arr * 1.60934;
    velo_rpm = [1500  9500];
    throttle_p = 1;
    ma_ideal = (D/C)*(rc/(rc-1)) * rho(1);
    D_cc = 1500;
    Vd = ((pi*bore^2*stroke)/4); %meter cubes
    V1  = Vd/(1-(1/rc));
    Qconst = combustion_eff * Qlhv/(p(1)*V1); 
    Ws = Cv*((T(3)-T(2)) - (T(4)-T(1)));
    power_needed = PowerReq(60) * 2
    velo1 = 60 * 1.61
    cylinders = 4
    for i = 1:length(N)
        if N(i) > 5000
            cylinders = 6
        end
        
        efficiency(i) = 1;
        RPM = [2100 15000];
        mech_eff = [ .9 .65];
        
        if N(i)<= 2100
            mechanical_eff = .9;
        else
            mechanical_eff = interp1(RPM, mech_eff, N(i), 'linear');
        end
        
        efficiency(i) = mechanical_eff;
         %Specific work per cylinder, J/kg
        %Ma(i) = Veff(i) * (D/C)*(rc/(rc-1)) * rho(1);
        power_loss_pump = (pump_work(i)/1000 * N(i)/120)/0.745699872   
        %throttle_p(i) = throttle(N(i), power_needed/mechanical_eff + power_loss_pump, Veff(i), ma_ideal, fuel_air, Qconst,  mechanical_eff);
        throttle_p(i) = 1;
        Ma(i) = min(Veff(i), throttle_p(i)) *  ma_ideal;%Mass of Air in each cylinder, kg
        Mf(i) = Ma(i)/fuel_air;
        Q(i) = Mf(i) * Qconst;
        %PL(i) = PumpingLoss(Ma(i))
        net_work(i) = (FiniteHeatWoschni(Q(i), N(i), Ma(i), 0)); %kJ   
        velo(i) = interp1(velo_rpm, velocity_kmhr, N(i), 'linear');
        %velo1 = 60 * 1.61
        %W_loss(i) = drag_power(velo(i)) * 120 /N(i)
        Wt(i) = efficiency(i) * net_work(i)  - pump_work(i)/1000;
        % divide by 120 sicne two cycles per ottocylce
        
        pl(i) =  1 - pump_work(i)/(1000 * net_work(i));
        
        Power(i) = cylinders * Wt(i) * N(i)/120; % watts
        P_rate(i) = .8 * Power(i); 
        
        imep(i) = Power(i)/1000*2*10^3/(D*1000*(N(i)/60));
        bmep(i) =  mechanical_eff*imep(i);
        %p_real(i) = ((P_total(i) - Pf(i))/(P_total(i))) * P_total(i)
        %P_hp(k,i) = P_total(i)/745.699872;
        Power_watts(i)=  Power(i);
        P_hp(i) = Power(i)/0.745699872;
        %P_hp(i) =  Power(i)/745.699872; %
        %velocity(i) = velocity_finder(P_hp(i));
        % SFC
        %SFC(i) = (C * Ma(i)/f)/(Wt(i)); %kg/kj
        SFC(i) = (Mf(i))/(Wt(i)); %kg/kj
        SFC_Converted(i) = SFC(i) * 3.6e6; %kg/Kw-hr
        Torque(i) = 1000 * Power(i)/((2*3.1415*N(i)/60));   
        fuel_eff(i) = 0;
        if N(i)>= velo_rpm(1) && N(i)<= velo_rpm(end)           
            % .9 because of transmission losses
            velo1 = velocity_finder(P_hp(i))
            fuel_eff(i) = .9 * (velo1/(SFC_Converted(i) * (Power_watts(i)/1000))) * 2.35215; % to get to MPG
        end
        
    end
%%
figure 

disp(N)
plot(N, P_hp)
hold on 
%plot(N,P_hp1)
title('Power with Cylinder Deactivation vs RPM')
ylabel('Power (hp)')
xlabel('RPM')
xlim([0 9300])
%legend('6 Cylinder', '4 Cylinder')
%xline(9300, 'LineWidth',3)

%%
figure 
hold on
plot(N, Torque(:,1))
title('Torque vs RPM')
xlabel('RPM')
ylabel('Torque, N*m')



%%
figure
plot(N, Veff)
title('Volumetric Efficiency vs RPM')
xlabel('RPM')
ylabel('Volumetric Efficiency');
xlim([0 9300])

%%
figure
plot(N, SFC_Converted)
title('SFC')
%%
figure
plot(N, fuel_eff)
ylabel('Fuel Efficiency (MPG)')
hold on
yyaxis right
plot(N,throttle_p)
title('Fuel Efficiency at 60 MPH, with Throttling')
xlabel('RPM')
ylabel('Throttle Needed to Cruise at 60 MPH')
legend('Fuel Efficiency', 'Throttle Percentage')
%%
figure
hold on
plot(N,throttle_p)
title('Throttle Needed to Cruise at 60 MPH')
xlabel('RPM')
ylabel('Throttle Percentage')

%%
figure
hold on
plot(N,velocity)
title('Velocity at Full Throttle (High Gear)')
xlabel('RPM')
ylabel('Velocity (MPH)')
%%
figure
plot(N, fuel_eff)
ylabel('Fuel Efficiency (MPG)')
hold on
yyaxis right
plot(N,velocity)
title('Fuel Efficiency & Velocity at Full Throttle (Highest Gear)')
xlabel('RPM')
ylabel('Velocity at Full Throttle (MPH)')
legend('Fuel Efficiency', 'Velocity')

%%
power_loss_pump = (pump_work.'/1000 .* N/120)/0.745699872  
figure
hold on
grid
plot(N,power_loss_pump)
title('Pump Power Loss (hp) vs RPM')
xlabel('RPM')
ylabel('Pump Power Loss (hp)')
xlim([0 9300])
ylim([min(power_loss_pump), 8])
