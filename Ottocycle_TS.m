%%
clear all
i = 1 
rc = 10
% Initial Conditions
p(1) = 101325;%pa
%p(1,i) = p(1)*1000; %pa
%T(1,i) = 288; %k
T(1) = 294; %k
gamma = 1.4;
R = 287; %j/kg, gas constant
%rc = 9; % compression ration
Qlhv = 42.7e6; %j/kg
Cp = R*gamma/(gamma-1); 
f = 14.7; % air to fuel ratio
D = .0015

% State 1
rho(1) = p(1)/(R*T(1)); 

% State 2
rho(2) = rho(1) * rc;
p(2) = p(1)* (rho(2)/rho(1))^gamma;
T(2) = p(2)/(rho(2)* R);

% State 3
T(3) = T(2) + .8*Qlhv/(Cp* f); %.8 is thermal efficeny of combustionn
rho(3) = rho(2);
p(3) = rho(3) * R* T(3);
T(3) = 3000

% State 4
rho(4) = rho(1);
p(4) = p(3)* (rho(4)/rho(3))^ (gamma-.1); %this is an estimate, could use refprop to find actual
T(4) = p(4)/(rho(4) * R);
j = 1;

% Work
%m^3
C = 6; % initial guess
Cv = Cp-R;
Ws = Cv*((T(3)-T(2)) - (T(4)-T(1)));
Ma = (D/C)*(rc/(rc-1)) * rho(1);
Wc = Ws*Ma;
Wt = Wc * C; %total work

%Power
N = 5000; % RPM, crusie
Power = Wt * N/120;
P_rate = .8 * Power; 
P_specific = Ws * N/120;
P_cylinder = P_specific * Ma;
P_total = P_cylinder * C;
P_hp = P_total/745.699872;
% SFC
SFC = (C * Ma/f)/Wt; %kg/kj
SFC_Converted = SFC * 3.6e9; %g/Kw-hr

 S(1)= refpropm('S','D',rho(1), 'P',p(1)/1e3, 'air.ppf');
 S(2)= refpropm('S','D',rho(2), 'P',p(2)/1e3, 'air.ppf');
 S(3)= refpropm('S','T',3000, 'P',p(3)/1e3, 'air.ppf');
 S(4)= refpropm('S','D',rho(4), 'P',p(4)/1e3, 'air.ppf');


%%
% Calculates non-linear S between states 2 and 3 
i = 2;
Tn =  T(2) + 10;
p_n(1) = Tn * R * rho(2)
S_n(1)= refpropm('S','D',rho(2), 'P',p_n(1)/1e3, 'air.ppf');
H_n(1)= refpropm('H','D',rho(2), 'P',p_n(1)/1e3, 'air.ppf');
for Tn  = Tn: 5 : T(3)
    T_n(i) = Tn;
    p_n(i) =   Tn * R * rho(2);
    S_n(i)= refpropm('S','D',rho(2), 'P',p_n(i)/1e3, 'air.ppf');
    H_n(i)= refpropm('H','D',rho(2), 'P',p_n(i)/1e3, 'air.ppf');
    i = i + 1; 
end
T_n(1) = T(2)
S_n(1) = S(2)

% Calculate non-linear S between 1 and 4
i = 2;
rho14 = rho(1)
T14 =  T(1) + 10;
p_14(1) = T14 * R * rho14
S_14(1)= refpropm('S','D',rho14, 'P',p_14(1)/1e3, 'air.ppf');
H_14(1)= refpropm('H','D',rho14, 'P',p_14(1)/1e3, 'air.ppf');
for T14  = T14: 5 : T(4)
    T_14(i) = T14;
    p_14(i) = T14 * R * rho14;
    S_14(i)= refpropm('S','D',rho14, 'P',p_14(i)/1e3, 'air.ppf');
    H_14(i)= refpropm('H','D',rho14, 'P',p_14(i)/1e3, 'air.ppf');
    i = i + 1; 
end
T_14(1) = T(1)
S_14(1) = S(1)

H(1)= refpropm('H','D',rho(1), 'P',p(1)/1e3, 'air.ppf');
H(2)= refpropm('H','D',rho(2), 'P',p(2)/1e3, 'air.ppf');
H(3)= refpropm('H','T',3000, 'P',p(3)/1e3, 'air.ppf');
H(4)= refpropm('H','D',rho(4), 'P',p(4)/1e3, 'air.ppf')

p_m = p/1e6
T(3) = T(2) + .8*Qlhv/(Cp* f);
scatter(S,T, 'MarkerEdgeColor','red')
hold on
scatter(S_n,T_n,'MarkerEdgeColor','red')
hold on
scatter(S_14,T_14, 'MarkerEdgeColor','red')
title('Ideal Otto Cycle TS Diagram')
xlabel('Entropy (kJ/kG-K)')
ylabel('Temperature (K)')