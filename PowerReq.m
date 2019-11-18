function Preq = PowerReq(vmph)
%Power Required at a certain Speed

Cr = 0.015;
Cd =  0.7;
M = 300; %kg
g = 9.81; %m/s^2
A = 1; % m^2
p = 1.23; %kg/m^3

F = Cr*M*g; %rolling Force
Q = 0.5*p*Cd*A; %Drag force
v2 = .44704*vmph; % m/s
Preq = (F+Q*v2^2)*v2/745.1; %hp



end