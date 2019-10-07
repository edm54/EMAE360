%Air intate
clear all
clc
close all

set(0, 'DefaultAxesFontWeight', 'normal', ...
      'DefaultAxesFontSize', 18, ...
      'DefaultAxesFontAngle', 'normal', ... 
      'DefaultAxesFontWeight', 'normal', ... 
      'DefaultAxesTitleFontWeight', 'bold', ...
      'DefaultAxesTitleFontSizeMultiplier', 1.2) ;
set(groot,'defaultLineLineWidth',3)
bore = 7.045/100;
r_bs =  1.1
stroke = bore/r_bs;
max_v_lift = .12 * bore; %meter
x = [-360 : 1: 360];
y1 = -1*(x.^2 -220*x - 1125);

max_y =max(y1);
% Recreating Heywood 225
intake_lift = -(max_v_lift/max_y)*(x.^2 -220*x - 1125);



hold on
y2 = -1.*(x.^2 + 215*x - 2250);
max_2 = max(y2);


exhuast_lift = -1 * (max_v_lift/max_2) .*(x.^2 + 215*x - 2250);
for i = 1:length(intake_lift)
    if intake_lift(i) < 0
        intake_lift(i) = 0;
    end
    
    if exhuast_lift(i) < 0
        exhuast_lift(i) = 0;
    end
end

plot(x, exhuast_lift*1000);
plot(x,intake_lift*1000);
title('Intake Valve Lift as a function of crank angle')
ylim([0,10]);
legend('Exhuast valve lift','Intake valve lift')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(gcf,'Valve Lift vs Crank.png')

for i = 1:length(intake_lift)
    if intake_lift(i) < 0
        intake_lift(i) = 0;
    end
    
    if exhuast_lift(i) < 0
        exhuast_lift(i) = 0;
    end
end

%%

% Fuerg 5.11
D_in = .333 * bore;
D_ex = .29 * bore ;
DS_in = .21 * D_in;
DS_ex = .24 * D_ex;

DV_in = 1.1 * D_in;
DV_ex = 1.1*D_ex;

DP_in = D_in;
DP_ex = D_ex;

% USe 221 nd 222 to find rest
w_ex = D_ex * (1.11 - 1)/2;
w_in = D_in * (1.11 - 1)/2; %P 221, 30 degrees

DM_in = DV_in - w_in
DM_ex = DV_ex - w_ex

B = 30 * pi/180;

for i = 1: length(intake_lift)
    if (w_in /(sin(B) * cos(B))) > intake_lift(i) && intake_lift(i)>0
        Am_in(i) = pi * intake_lift(i)*cos(B) * (DV_in - 2*w_in + (intake_lift(i)/2)*sin(2* B));
        disp('I 1')
    elseif (((((DP_in^2 - DS_in^2)/(4*(DM_in)))^2 - w_in^2)^(1/2) + w_in * tan(B)) >= intake_lift(i) ...
            && (w_in /(sin(B) * cos(B))) < intake_lift(i) ...
            && intake_lift(i)>0)
        Am_in(i) = pi *  DM_in * sqrt((intake_lift(i) - w_in * tan(B))^2 + w_in^2);
        disp('I2')
    elseif intake_lift(i)>0
        Am_in(i) = pi * (DP_in^2 - DS_in^2)/4;
        disp('I3')
    else
        Am_in(i) = 0;

    end
end

for i = 1: length(exhuast_lift)
    if (w_ex /(sin(B) * cos(B))) > exhuast_lift(i) && exhuast_lift(i)>0
        Am_ex(i) = pi * exhuast_lift(i)*cos(B) * (DV_ex - 2*w_ex + (exhuast_lift(i)/2)*sin(2* B));

    elseif sqrt(((DP_ex^2 - DS_ex^2)/(4*(DM_ex)))^2  - w_ex^2) + w_ex * tan(B) >= exhuast_lift(i) && exhuast_lift(i)>0
        Am_ex(i) = pi *  (DM_ex) * sqrt((exhuast_lift(i) - w_ex * tan(B))^2 + w_ex^2);     

    elseif exhuast_lift(i)>0
        Am_ex(i) = pi * (DP_ex^2 - DS_ex^2)/4;

    else
         Am_ex(i) = 0;
    end
    
end
 
figure
hold on 
plot(x,Am_in)
plot(x, Am_ex)

legend('Intake', 'Exhuast')
title('Flow Area (m^2) as a Function of Crank Angle')
legend('Intake Flow Area', 'Exhuast Flow Area')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(gcf,'Flow Area vs Crank.png')


%% ex
cd_ex = .5 + .5 .* [  1.2 1.5 1.6 1 .4 .2 0]./3.2


lv_dv_ex = [.1 .15 .2 .25 .3 .35 .4]
plot(lv_dv_ex, cd_ex)
fit = polyfit(lv_dv_ex, cd_ex,  3)
val = polyval(fit, .1:.01:.3)
hold on
plot(.1:.01:.3, val)
ylim([0 1])
xlim([0 .3])

%% intake Cd vs lv/dv
figure

l_in1 = [.05 .08 .1]
c_1 = [.58 .62 .67]

l_in2 = [.1 .11]
c_2 = [.67 .55] 

l_in3 = [.11 .14 .18]
c_3 = [.55 .59 .68]

l_in4 = [.18 .21 .24 .28]
c_4 = [.68 .55 .53 .48]

fit1 = polyfit(l_in1, c_1,  3)
val1 = polyval(fit1, .05:.01:.1)
plot(.05:.01:.1, val1)
hold on
fit2 = polyfit(l_in2, c_2,  1)
val2 = polyval(fit2, l_in2)
plot(l_in2, val2)

fit3 = polyfit(l_in3, c_3,  2)
val3 = polyval(fit3, .11:.01:.18)
plot(.11:.01:.18, val3)

fit4 = polyfit(l_in4, c_4,  1)
val4 = polyval(fit4, .18:.01:.35)
plot(.18:.01:.35, val4)

title('Intake Cd vs LV/DV')
%%


for i = 1: length(intake_lift)
    lv_dv(i) = intake_lift(i)/DV_in;
    if lv_dv(i) < .1 && lv_dv(i) > 0
        CD_in(i) = polyval(fit1, lv_dv(i));
    elseif lv_dv(i) <= .11 && lv_dv(i) > 0
        CD_in(i) = polyval(fit2, lv_dv(i));
    elseif lv_dv(i) <= .18 && lv_dv(i) > 0
        CD_in(i) = polyval(fit3, lv_dv(i));
    elseif lv_dv(i) > 0 
        CD_in(i) = polyval(fit4, lv_dv(i));
    else
        CD_in(i) = 0;
    end
    
    
    lv_dv2(i) = exhuast_lift(i)/DV_ex;
    CD_ex(i) = polyval(fit, lv_dv2(i));
    if exhuast_lift(i) == 0
        CD_ex(i) = 0;
    end
    
    
end 
figure
plot(x, CD_in)
hold on
plot(x, CD_ex)
title('CD as a Function of Crank Angle')
xlabel('Crank angle')
ylabel("Discharge Coeff")
legend('Intake', 'Exhuast')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(gcf,'CD vs Crank.png')
%%
N = 5000
N = 15000
for theta = 0:360
    time = (theta/360) * (1/(N/60));
end

%% Exhuast
N = 5000
To = 1.9286e+03 ;%k
gamma = calc_gamma(To);
P_external_e = 101325; %pa

R = 287;
Co = sqrt(gamma * R * To);
P_cylinder_e = 6.6417e+05; % Pascals
rho_cylinder = P_cylinder_e/(R * To); %Kg/m3
Ma = 3.335683867042752e-04; %Kg
% D_ex is diamter
%x = [-360 : 1: 360];
theta = [-360 : 360];
time = (theta./360) .* (1/(N/60)); %seconds
for i = 1:length(theta)
    v(i) = calc_volume(theta(i));
end

temp = To;
pressure = P_cylinder_e * ones(1,length(theta));
mass = Ma;
d_mass_d_time = 0;
delta_t = time(2) - time(1);


for i = 1:length(theta)-1
    lift = exhuast_lift(i);
    min_area = Am_ex(i) ; %meters
    gamma = calc_gamma(temp(i));
    Cd = CD_ex(i);
    if pressure(i) > P_external_e 
        m_dot(i) = -1*Cd * rho_cylinder * min_area * Co * sqrt((2/(gamma-1)) * ((P_external_e/P_cylinder_e)^(2/gamma) - ...
        (P_external_e/P_cylinder_e)^((gamma+1)/gamma))) ;
    else
        % Flow into cylinder from exhuast
        m_dot(i) =   Cd * rho_cylinder * min_area * Co * sqrt((2/(gamma-1)) * ((P_external_e/P_cylinder_e)^(2/gamma) - ...
        (P_external_e/P_cylinder_e)^((gamma+1)/gamma)));
    end
    
    % Choked flow?
    rho(i) = pressure(i)/(R * temp(i));
    velocity(i) = abs(m_dot(i))/(rho(i) * min_area);
    sound_speed(i) = sqrt(gamma * R *  temp(i));
        
    if velocity(i) > sound_speed(i)
        m_dot(i) = Cd * rho_cylinder * min_area * Co *(2/(gamma + 1)) ^ ((gamma + 1)/(2 *(gamma - 1)));
        disp('hoke')
        disp(i)
    end
    d_mass_d_time = m_dot(i);
    mass(i + 1) = mass(i) + delta_t * d_mass_d_time;
    temp(i+1) = temp(i);
    if m_dot(i) ~= 0
        pressure(i + 1) =  mass(i + 1) * R * temp(i+1)/v(i+1);
    end
    

end

%U = sqrt(2* gamma/(gamma - 1) * (Po/rho0) * (1 - (Pv/Po)^((gamma - 1)/gamma)))

%plot( theta(1 :720), rho)
%plot( theta, mass)
figure
plot( theta, pressure)
figure
plot(theta, mass)
