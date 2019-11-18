function [P_needed] = drag_power(kmhr)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
cd = .7;
rho = 1.2; %kg/m3
A = 1; %m2

ms_velo = 0.277778 * kmhr;
velo = ms_velo;% m/s
D = .5 * rho * A * cd * velo^2 ;
roll_resistance = .015 * 181 * 9.81;
F = D + roll_resistance;

P_needed = F * velo;
%P_needed_hp = P_needed/746 
end

