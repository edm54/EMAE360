N = [800, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 9300];
Vtdc = 2.781e-5; %m^3
V1 = Vtdc*10;
Qlhv = 44e6; %j/kg from Heywood App D
p1 = 101325;

Qin = combustionEff(1)*Qlhv/(p1*V1); 
rho = zeros(length(N),1);
temp = zeros(length(N),1);
nv = volumetric_efficiency(N);
for i = 1:length(N)
    ma = 3.335e-4; %kg
    ma = ma*nv(i);
    Qin = ma*combustionEff(1)*Qlhv/(p1*V1); 
    rho(i) = ma/Vtdc;
    temp(i) = Tcombust(Qin,N(i),ma,0);
end

N = 5000;
ma = 0.5 * 3.335e-4; %kg
Qin = ma*combustionEff(1)*Qlhv/(p1*V1); 
rho=ma/Vtdc;
temp = Tcombust(Qin,N,ma,0);