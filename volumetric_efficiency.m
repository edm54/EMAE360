function [V_eff] = volumetric_efficiency(N)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% all
%close all
%% Prepare Pllots
set(0, 'DefaultAxesFontWeight', 'normal', ...
      'DefaultAxesFontSize', 18, ...
      'DefaultAxesFontAngle', 'normal', ... 
      'DefaultAxesFontWeight', 'normal', ... 
      'DefaultAxesTitleFontWeight', 'bold', ...
      'DefaultAxesTitleFontSizeMultiplier', 1.2) ;
set(groot,'defaultLineLineWidth',3)

bore = 7.045/100;
r_bs =  1.1;
stroke = bore/r_bs;
%max_in_lift = .12 * bore; %meter
%max_ex_lift = .12 * bore;
max_in_lift = .007; %meter
max_ex_lift = .008459; % meter

theta = [-360 : .1: 360];
y1 = -1*(theta.^2 -220*theta - 1125);

open_in = -5;
in_dur = 220;

close_in = in_dur + open_in;
intake_lift = -1 * (theta - open_in).* (theta - close_in); 

open_ex = -220;
ex_dur = 230;
close_ex = ex_dur + open_ex;

exhuast_lift = -1 * (theta - open_ex).* (theta - close_ex); 
y_out_max = max(exhuast_lift)
y_in_max = max(intake_lift)

intake_lift = (max_in_lift/y_in_max).* intake_lift
exhuast_lift = (max_ex_lift/y_out_max) .* exhuast_lift

for i = 1:length(intake_lift)
    if intake_lift(i) < 0
        intake_lift(i) = 0;
    end
    if exhuast_lift(i) < 0
        exhuast_lift(i) = 0;
    end
end

plot(theta, exhuast_lift*1000);
hold on
plot(theta,intake_lift*1000);
title('Intake Valve Lift as a function of crank angle')
xlabel('Crank Angle')
ylabel('Lift (mm)')
ylim([0,10]);
legend('Exhuast valve lift','Intake valve lift')

%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%saveas(gcf,'Valve Lift vs Crank.png')
    
% Fuerg 5.11
% 2 times 
 %D_in = .33 * bore;
 %D_ex = .29 * bore ;

 D_in = .026;
 D_ex = .0225 ;
%D_ex = .019 ;

DS_in = .21 * D_in;
DS_ex = .24 * D_ex;

DV_in = 1.1 * D_in;
DV_ex = 1.1*D_ex;

DP_in = D_in;
DP_ex = D_ex;

% USe 221 nd 222 to find rest
w_ex = D_ex * (1.11 - 1)/2;
w_in = D_in * (1.11 - 1)/2; %P 221, 30 degrees

DM_in = DV_in - w_in;
DM_ex = DV_ex - w_ex;

B = 30 * pi/180;
%% Flow Area
for i = 1: length(intake_lift)
    if (w_in /(sin(B) * cos(B))) > intake_lift(i) && intake_lift(i)>0
        Am_in(i) = pi * intake_lift(i)*cos(B) * (DV_in - 2*w_in + (intake_lift(i)/2)*sin(2* B));
       %disp('I 1')
    elseif (((((DP_in^2 - DS_in^2)/(4*(DM_in)))^2 - w_in^2)^(1/2) + w_in * tan(B)) >= intake_lift(i) ...
            && (w_in /(sin(B) * cos(B))) < intake_lift(i) ...
            && intake_lift(i)>0)
        Am_in(i) = pi *  DM_in * sqrt((intake_lift(i) - w_in * tan(B))^2 + w_in^2);
        %disp('I2')
    elseif intake_lift(i)>0
        Am_in(i) = pi * (DP_in^2 - DS_in^2)/4;
        %disp('I3')
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
% Or flat at .8 percent 
figure
hold on 

plot(theta,Am_in)
plot(theta, Am_ex)

legend('Intake', 'Exhuast')
title('Flow Area (m^2) as a Function of Crank Angle')
legend('Intake Flow Area', 'Exhuast Flow Area')
ylabel('Flow Area (M^2)')
xlabel('Crank Angle')
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%saveas(gcf,'Flow Area vs Crank.png')
figure

%% ex
cd_ex = .5 + .5 .* [  1.2 1.5 1.6 1 .4 .2 0]./3.2;


lv_dv_ex = [.1 .15 .2 .25 .3 .35 .4];
plot(lv_dv_ex, cd_ex);
fit = polyfit(lv_dv_ex, cd_ex,  3);
val = polyval(fit, .1:.01:.3);
hold on
plot(.1:.01:.3, val);
ylim([0 1]);
xlim([0 .3]);

%% intake Cd vs lv/dv
figure

l_in1 = [.05 .08 .1];
c_1 = [.58 .62 .67];

l_in2 = [.1 .11];
c_2 = [.67 .55] ;

l_in3 = [.11 .14 .18];
c_3 = [.55 .59 .68];

l_in4 = [.18 .21 .24 .28];
c_4 = [.68 .55 .53 .48];

fit1 = polyfit(l_in1, c_1,  2);
val1 = polyval(fit1, .05:.01:.1);
plot(.05:.01:.1, val1);
hold on
fit2 = polyfit(l_in2, c_2,  1);
val2 = polyval(fit2, l_in2);
plot(l_in2, val2);

fit3 = polyfit(l_in3, c_3,  2);
val3 = polyval(fit3, .11:.01:.18);
plot(.11:.01:.18, val3);

fit4 = polyfit(l_in4, c_4,  1);
val4 = polyval(fit4, .18:.01:.35);
plot(.18:.01:.35, val4);

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
plot(theta, CD_in)
hold on
plot(theta, CD_ex)
title('CD as a Function of Crank Angle')
xlabel('Crank angle')
ylabel("Discharge Coeff")
legend('Intake', 'Exhuast')
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%saveas(gcf,'CD vs Crank.png')

%% Exhuast
ind = 1;

%N = 800:25:9400;


for j = 1 : length(N)
    V_eff(ind) = 1;
    To = 300 ;%k
    Tf = 1.9286e+03 ;%k
    gamma = calc_gamma(To);
    P_external_e = 100000;
    P_external_in = 101325;
    R = 287;
    Co = sqrt(gamma * R * To);
    P_cylinder_e = 6.6417e+05; 
    P_cylinder_e = 4.6417e+05;% Pascals
    rho_cylinder = P_cylinder_e/(R * To); %Kg/m3
    Ma = 3.335683867042752e-04; %Kg
    time = (theta./360) .* (1/(N(j)/60)); %seconds
    for i = 1:length(theta)
        v(i) = calc_volume(theta(i));
    end

    temp = Tf;
    pressure = P_cylinder_e * ones(1,length(theta));
    mass = Ma;
    d_mass_d_time = 0;
    delta_t = time(2) - time(1);
    rho = rho_cylinder;
    m1 = Ma;

    for i = 1:length(theta)-1
        rho(i) = pressure(i)/(R * temp(i));
        m1(i) = rho(i)*v(i);
        ex_lift = exhuast_lift(i);
        in_lift = intake_lift(i);
        min_area_ex = Am_ex(i); %meters
        min_area_in = Am_in(i);
        gamma = calc_gamma(temp(i));
        Cd_ex = CD_ex(i);
        Cd_in = CD_in(i);
        m_dot(i) = 0;
        ex_close = 10;
        in_open = -5;
        if theta(i) < ex_close
            r1(i) = P_external_e/pressure(i);
            %  Check flow into cylinder
            if pressure(i) > P_external_e  
                if (min_area_ex > 0)
                    m_dot(i) = m_dot(i) + -2 * (((Cd_ex * min_area_ex * pressure(i))/sqrt(R * temp(i)))...
                        * r1(i)^(1/gamma) * sqrt(((2 * gamma)/(gamma - 1)) * (1 - r1(i)^((gamma -1)/gamma))));
                else
                    m_dot(i) = 0;
                end
                
                % Check choked flow
                if r1(i) <= (2/(gamma + 1))^(gamma/(gamma - 1))/15 && m_dot(i)<0                   
                     m_dot(i) = -2 * Cd_ex * pressure(i) * min_area_ex * (gamma /Co) *(2/(gamma + 1))...
                        ^ ((gamma + 1)/(2 *(gamma - 1)));  
                end
            else
                
                if (min_area_ex > 0)
                % Flow into cylinder from exhuast
                    m_dot(i) = 2 * Cd_ex * rho(i) * min_area_ex * Co * sqrt((2/(gamma-1))...
                        * ((P_external_e/pressure(i))^(-2/gamma) - ...
                    (r1(i))^(-1*(gamma+1)/gamma)));
                else
                    m_dot(i) = 0;
                end
            end
            
        end            
        if theta(i) >= in_open 
               
            % Track Volumetric Efficicency
            if theta(i) == -5
                 mass_start = mass(i);
            end
            
            r2(i) = pressure(i)/P_external_in;
            % Flow into the cylinder
            if pressure(i) < P_external_in
                m_dot(i) = m_dot(i) + 2 * (((Cd_in * min_area_in * P_external_in)/sqrt(R * temp(i)))...
                    * r2(i)^(1/gamma) * sqrt(((2 * gamma)/(gamma - 1)) * (1 - r2(i)^((gamma -1)/gamma))));
                    
                % Check choked flow
                if r2(i) <= ((2/(gamma + 1))^(gamma/(gamma - 1)))/15            
                     m_dot(i) =  m_dot(i) + 2  * Cd_in * rho_cylinder * min_area_in * (gamma /Co) *(2/(gamma + 1))...
                        ^ ((gamma + 1)/(2 *(gamma - 1)));
                end
            else
                % Flow out of cylinder into inlet
                    m_dot(i) = m_dot(i) + -2 * Cd_in * rho(i) * min_area_in * Co * sqrt((2/(gamma-1))...
                        * ((r2(i))^(-2/gamma) - (r2(i))^(-1*(gamma+1)/gamma))); 
            end
        end   
        d_mass_d_time = m_dot(i);
        % Update mass in cylinder
        mass(i + 1) = max(mass(i) + delta_t * d_mass_d_time, .000001);
        if mass(i+1) < 0
            mass(i+1) = .00001;
        end

       if theta(i) <= -5
            temp(i+1) = temp(i);
       else
            temp(i+1) = (mass(i) * temp(i) + To * (mass(i+1) - mass(i)))/mass(i+1);
       end
       if abs(m_dot(i)) > .00001
            pressure(i + 1) =  mass(i + 1) * R * temp(i+1)/v(i+1);
       end

    end
    % Calculate individual Volumetric Efficiencies
    veff2(ind) =  mass(end)/mass(1);
    veff1(ind) = (1 - mass_start/mass(1));
    V_eff(ind) = V_eff(ind) * veff1(ind) * veff2(ind);
    
    ind = ind + 1;
    %m_dot * temp * dt + temp in * mass
end

    figure
    hold on
    plot(N, V_eff)
    plot(N, veff1)
    plot(N, veff2)
    grid 
    title('Volumetric Efficency vs RPM')
    legend('Total', 'Exhuast', 'Intake')
    
    
    
   disp(max(V_eff))
   disp(max(veff1))
   disp(max(veff2))
    
    %%
    
    figure
    plot(theta, mass)
    title('Air mass in cylinder during exhuast as a function of crank angle')
    xlabel('Crank Angle')
    ylabel('Mass of Air (kg)')
    
    figure
    plot(theta, temp)
    title('temp')
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    %saveas(gcf,'Air mass exhaust vs crank.png')

    figure
    plot(theta(1:end-1), m_dot)
    title('M dot during exhuast as a function of crank angle')
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    %saveas(gcf,'M dot vs Crank Angle.png')

    
 %%   
%     %U = sqrt(2* gamma/(gamma - 1) * (Po/rho0) * (1 - (Pv/Po)^((gamma - 1)/gamma)))
% 
%     %plot( theta(1 :720), rho)
%     %plot( theta, mass)
%     figure
%     plot(theta(1:720), sound_speed)
%     hold on
%     plot(theta(1:720), velocity)
%     title('sound speed & velo vs theta')
%%
%     figure
%     plot(theta(1:end-1), rho)
%     title('rho vs theta')
% 
%   
% %%
%     figure
%     plot( theta, pressure)
%     title('Pressure of cylinder during exhuast as a function of crank angle')
%     line([theta(1),theta(end)],[101000,101000])
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    %saveas(gcf,'P exhuast vs crank.png')
%%

%     figure
%     plot(theta, mass)
%     title('Air mass in cylinder during exhuast as a function of crank angle')
%     %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%     %saveas(gcf,'Air mass exhaust vs crank.png')
% 
%     figure
%     plot(theta(1:end-1), m_dot)
%     title('M dot during exhuast as a function of crank angle')
%     %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%     %saveas(gcf,'M dot vs Crank Angle.png')


 
 %%
 %N = 800: 25: 10000
figure
hold on
plot(N, V_eff)
plot(N, veff1)
plot(N, veff2)
grid 
title('Volumetric Efficency vs RPM')
legend('Exhuast', 'Intake', 'Total')
    %}
   