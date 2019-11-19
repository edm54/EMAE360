af = 14.9;
Vtdc = 2.781e-5; %m^3
V1 = Vtdc*10;
Qlhv = 44e6; %j/kg from Heywood App D
p1 = 101325;
X_1CO = [.02389,.023868,.0238, .023579,.023733,.023715,.023709,.023726,.023777,.023867,.023902];
X_1NO = [.0081928, .0081908,.0081845, .0081808,.0081783,.0081767,.0081761,.0081777,.0081824,.0081907,.0081939];
N1 = [800,1000,2000,3000,4000,5000,6000,7000,8000,9000,9300];
CO = polyfit(N1,X_1CO,1);
NO = polyfit(N1,X_1NO,1);

N = 800:100:9300;
co_mpkm = zeros(length(N),1);
no_mpkm = zeros(length(N),1);
xco = zeros(length(N),1);
xno = zeros(length(N),1);
vel = zeros(length(N),1);
RPM = [2100 15000];
mech_eff = [ .9 .65];  
for i = 1:length(N)
    ma = 3.335e-4; %kg
    [nv,pump_work] = volumetric_efficiency(N(i));
    close all
    ma = ma*nv;
    mf = ma/14.9;
    Qin = mf*combustionEff(1)*Qlhv/(p1*V1); 
    net_work = 6 * (FiniteHeatWoschni(Qin, N(i), ma, 0)); %kJ   
    if N(i)<= 2100
        mechanical_eff = .9;
    else
        mechanical_eff = interp1(RPM, mech_eff, N(i), 'linear');
    end
    efficiency = mechanical_eff;
    Wt(i) = efficiency * net_work/6  - pump_work/1000;

    Power(i) = 6 * Wt(i) * N(i)/120; % watts
    P_hp(i) = Power(i)/0.745699872;
    vel(i) = velocity_finder(P_hp(i)); %mph
end
for i = 1:length(N)
    nv = volumetric_efficiency(N(i));
    close all
    ma = 0.0003335; %kg, per cylinder
    ma = ma*nv;
    af = 14.9;
    mt = ma+ma/af; %kg
    mt = mt * 1000; %g
    v = vel(i)/2.237; %m/s
    mpkm = mt*(1/2)*(N(i)/60)*(1/v)*1000; %total mass of exhaust per km
    xco(i) = polyval(CO,N(i));
    xno(i) = polyval(NO,N(i));
    co_mpkm(i) = xco(i)*mpkm*6;%times 6 to account for 6 cylinders
    no_mpkm(i) = xno(i)*mpkm*6;
end

ReqPerCO = 1 - 12 / max(co_mpkm);
ReqPerNO = 1 - 0.8 / max(no_mpkm);

plot(N,co_mpkm,N,no_mpkm)
title(["Emissions vs RPM"; "Stoichiometric"; "Max Velocity (WOT)"])
xlabel("RPM")
ylabel("g of gas per km")
legend("CO (g/km)", "NO+HC (g/km)", "Location", "northwest")