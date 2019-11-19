function [throttle1] = throttle(RPM,power_needed, Veff, Ma_opt, fuel_air, Q_const, m_eff, cylinders)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
P_hp1  = 160;
throttle1 = 1;
while abs(P_hp1 - power_needed) > 1
    Ma1 = min(Veff, throttle1) * Ma_opt; %Mass of Air in each cylinder, kg
    Mf1 = Ma1/fuel_air;
    Q1 = Mf1 * Q_const;
    [v, pump_actual] = volumetric_efficiency_p(RPM, min(101000, throttle1*101000));
    net_work1 = FiniteHeatWoschni(Q1, RPM, Ma1, 0); %kJ
    work_tot = m_eff * net_work1 - pump_actual/1000;
    Power1 = cylinders * work_tot * RPM/120;
    P_hp1 = Power1/0.745699872
    throttle1 = throttle1 - .01 * (P_hp1 - power_needed);
    
end


end

