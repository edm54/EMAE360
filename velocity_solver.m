%Power Required at a certain Speed
clear
Cr = 0.015;
Cd = 0.4;
M = 300; %kg
g = 9.81; %m/s^2
A = 1; % m^2
p = 1.23; %kg/m^3
iter = 1;

F = Cr*M*g;
Q = 0.5*p*Cd*A;
vmph = 70; %mph, initial guess
v2 = .44704*vmph; % m/s
v_conv = 1;
Pa = 55; %hp

while v_conv >0.001
vi = v2;
Preq = (F+Q*vi^2)*vi/745.1; %hp
Pi = Preq - Pa;
P_prime = F+3*Q*(vi^2);

v2 = vi - Pi/P_prime;
v_conv = v2 - vi;
iter = iter+1;
end
v2 = v2*2.24; %mph

disp('mph')
disp(v2)

