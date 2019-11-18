%https://canadamotoguide.com/2016/05/04/motorcycle-aerodynamics/
cd = .7;
rho = 1.2; %kg/m3
A = 1; %m2
kmhr = 240;
ms_velo = 0.277778 * kmhr;
velo = ms_velo;% m/s
D = .5 * rho * A * cd * velo^2 
roll_resistance = .015 * 181 * 9.81
F = D + roll_resistance

P_needed = F * velo
P_needed_hp = P_needed/746 

p_hp = 180;

power_watts = p_hp * 746;
