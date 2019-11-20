%% Calculating Emissions at a certain velocity and throttle
clear
Vtdc = 2.781e-5; %m^3
V1 = Vtdc*10;
Qlhv = 44e6; %j/kg from Heywood App D
p1 = 101325;
%{
X_1CO = [.030263,.030241,.030171,.03013,.030103,.030085,.030079,.030097,.030148,.03024,.030275];
X_1NO = [.0077692,.0077663,.0077572,.0077518,.0077483,.0077459,.0077451,.0077474,.0077542,.0077662,.0077708];
N1 = [800,1000,2000,3000,4000,5000,6000,7000,8000,9000,9300];
CO = polyfit(N1,X_1CO,1);
NO = polyfit(N1,X_1NO,1);
N = 800:100:9300;
co1 = polyval(CO,N);
no1 = polyval(NO,N);



hold on
plot(N1,X_1NO,'o')
plot(N1,X_1CO, 'o')
plot(N,no1,N,co1)
xlabel("RPM")
ylabel("Mass fraction of total exhaust mass")
title("Mass fraction of CO and NOx vs RPM")
legend("NOx, points", "CO, points", "NOx, curve fit", "CO, curve fit")
hold off
%}
N = 3000;
throttle = 0.525;
vel = 60; %mph

ma = 0.0003335; %kg, per cylinder
ma = ma*throttle;
af = 14.7;
mt = ma+ma/af; %kg
mt = mt * 1000; %g
v = vel/2.237; %m/s
mpkm = mt*(1/2)*(N/60)*(1/v)*1000; %total grams of exhaust per km
%xco = polyval(CO,N); %mass fraction at WOT
%xno = polyval(NO,N); %Mass fraction at WOT
%xco = 0.01328; %mass fraction at a throttle of 0.355
%xno = 0.00267; %mass fraction at a throttle of 0.355
xco = 0.0127;
xno = 0.00236;
co_mpkm = xco*mpkm*4;%times 6 to account for 6 cylinders
no_mpkm = xno*mpkm*4;


ReqPerNO60 = 1 - 0.8 / max(no_mpkm);


