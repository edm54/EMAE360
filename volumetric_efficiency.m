function [V_eff] = volumetric_efficiency(N)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

bore = 7.045/100;
r_bs =  1.1;
stroke = bore/r_bs;
max_in_lift = .007 %meter
max_ex_lift = .008459 % meter

theta = [-360 : .25: 360];
y1 = -1*(theta.^2 -220*theta - 1125);

open_in = -5
in_dur = 220
close_in = in_dur + open_in
intake_lift = -1 * (theta - open_in).* (theta - close_in); 


open_ex = -220
ex_dur = 230
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


%%
%(theta - open) * (theta - close)
%find max and normalize it
%{
max_y =max(y1);
% Recreating Plot from Heywood 225
intake_lift = -(max_in_lift/max_y)*(theta.^2 -220*theta - 1125);

hold on
y2 = -1.*(theta.^2 + 215*theta - 2250);
max_2 = max(y2);


exhuast_lift = -1 * (max_ex_lift/max_2) .*(theta.^2 + 215*theta - 2250);
for i = 1:length(intake_lift)
    if intake_lift(i) < 0
        intake_lift(i) = 0;
    end
    if exhuast_lift(i) < 0
        exhuast_lift(i) = 0;
    end
end
%}

%%
% Look at 5000 hieght and AOC 
% Max areas:
   % In: 2.6 cm
   % Ex: 2.25 cm
   % Min is 1 cm 
   % .8459 cm of lift
    % floor of .1 cm
    
% Fuerg 5.11
% 2 times 
% D_in = .37 * bore;
% D_ex = .32 * bore ;

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

%% ex
cd_ex = .5 + .5 .* [  1.2 1.5 1.6 1 .4 .2 0]./3.2;

lv_dv_ex = [.1 .15 .2 .25 .3 .35 .4];
fit = polyfit(lv_dv_ex, cd_ex,  3);
val = polyval(fit, .1:.01:.3);

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
fit2 = polyfit(l_in2, c_2,  1);
val2 = polyval(fit2, l_in2);
fit3 = polyfit(l_in3, c_3,  2);
val3 = polyval(fit3, .11:.01:.18);
fit4 = polyfit(l_in4, c_4,  1);
val4 = polyval(fit4, .18:.01:.35);
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
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%saveas(gcf,'CD vs Crank.png')

%% Exhuast

%N = 800:25:9400;
%N = 800:100:15000;
ind = 1
for j = 1 : length(N)
    V_eff(ind) = 1;
    %choke(ind) = 0
    %reverse(ind) = 0
    To = 300 ;%k
    Tf = 1.9286e+03 ;%k
    gamma = calc_gamma(To);
    P_external_e = 100000;
    %P_external_e = 50000;
    %P_external_e = 101325
    % 101325; %pa TODO TRY
    P_external_in = 101325;
    %P_external_in = 651325;
    R = 287;
    Co = sqrt(gamma * R * To);
    P_cylinder_e = 6.6417e+05; 
    P_cylinder_e = 4.6417e+05;% Pascals
    %P_cylinder_e = 1.5e+05; % Pascals
    rho_cylinder = P_cylinder_e/(R * To); %Kg/m3
    Ma = 3.335683867042752e-04; %Kg
    % D_ex is diamter
    %x = [-360 : 1: 360];
    theta = [-360 : .25 :  360];
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
            if pressure(i) > P_external_e 
                 
                %could make rho consrtant i gues?
                if (min_area_ex > 0)
                    %m_dot(i) = -1*Cd * rho(i) * min_area * Co * sqrt((2/(gamma-1)) * ((P_external_e/pressure(i))^(2/gamma) - ...
                    %(P_external_e/pressure(i))^((gamma+1)/gamma))) ;
                    m_dot(i) = m_dot(i) + -2 * (((Cd_ex * min_area_ex * pressure(i))/sqrt(R * temp(i)))...
                        * r1(i)^(1/gamma) * sqrt(((2 * gamma)/(gamma - 1)) * (1 - r1(i)^((gamma -1)/gamma))));
                else
                    m_dot(i) = 0;
                end

                % Choked flow?
                %velocity(i) = abs(m_dot(i))/(rho(i) * min_area_ex);
                %sound_speed(i) = sqrt(gamma * R *  temp(i));

                %if velocity(i) > sound_speed(i)
                if r1(i) <= (2/(gamma + 1))^(gamma/(gamma - 1))/15 && m_dot(i)<0
                    %m_dot(i) = -1 * Cd * rho_cylinder * min_area * Co *(2/(gamma + 1))...
                       % ^ ((gamma + 1)/(2 *(gamma - 1)));
                     m_dot(i) = -2 * Cd_ex * pressure(i) * min_area_ex * (gamma /Co) *(2/(gamma + 1))...
                        ^ ((gamma + 1)/(2 *(gamma - 1)));

                   %m_dot(i) = -1 * ((CD_ex(i) * min_area * P_external_e)/sqrt(R * temp(i)))...
                   % *sqrt(gamma) * (2/(gamma + 1))^((gamma + 1)/(2*(gamma - 1))) ;
                    %%disp('choke')
                    %choke(ind) = choke(ind) + 1;
                    %disp(N(j))
                end
            else
                
                if (min_area_ex > 0)
                % Flow into cylinder from exhuast
                    m_dot(i) = 2 * Cd_ex * rho(i) * min_area_ex * Co * sqrt((2/(gamma-1))...
                        * ((P_external_e/pressure(i))^(-2/gamma) - ...
                    (r1(i))^(-1*(gamma+1)/gamma)));
                    %reverse(ind) = reverse(ind) + 1;
                    %disp('reverse')
                else
                    m_dot(i) = 0;
                end
            end
            
        end
        
%             
%         elseif 
%             %overlap
%             % overlaps -5 and 15
%             same as above
%             may not needif statement   
            
        if theta(i) >= in_open 
            %m_dot(i) = 0;
            if theta(i) == -5
                 mass_start = mass(i);
            end
            
            r2(i) = pressure(i)/P_external_in;
            % flow in
            %if min_area_in <= 0
                 %m_dot(i) = m_dot(i);
        %else
            if pressure(i) < P_external_in
                 
                m_dot(i) = m_dot(i) + 2 * (((Cd_in * min_area_in * P_external_in)/sqrt(R * temp(i)))...
                    * r2(i)^(1/gamma) * sqrt(((2 * gamma)/(gamma - 1)) * (1 - r2(i)^((gamma -1)/gamma))));
                    
                    %disp('not reverse ')
                    %disp(theta(i))
                % Choked flow?
                %velocity(i) = abs(m_dot(i))/(rho(i) * min_area_ex);
                %sound_speed(i) = sqrt(gamma * R *  temp(i));

                %if velocity(i) > sound_speed(i)
                if r2(i) <= ((2/(gamma + 1))^(gamma/(gamma - 1)))/15 
                    %m_dot(i) = -1 * Cd * rho_cylinder * min_area * Co *(2/(gamma + 1))...
                       % ^ ((gamma + 1)/(2 *(gamma - 1)));
                     m_dot(i) =  m_dot(i) + 2  * Cd_in * rho_cylinder * min_area_in * (gamma /Co) *(2/(gamma + 1))...
                        ^ ((gamma + 1)/(2 *(gamma - 1)));
                    
                   choke(ind) = choke(ind) + 1;
                   %m_dot(i) = -1 * ((CD_ex(i) * min_area * P_external_e)/sqrt(R * temp(i)))...
                   % *sqrt(gamma) * (2/(gamma + 1))^((gamma + 1)/(2*(gamma - 1))) ;
                    %disp('choke')
                    %disp(j)
                end
            else
                % Flow out of cylinder into inlet
                    m_dot(i) = m_dot(i) + -2 * Cd_in * rho(i) * min_area_in * Co * sqrt((2/(gamma-1))...
                        * ((r2(i))^(-2/gamma) - (r2(i))^(-1*(gamma+1)/gamma)));
                    %reverse(ind) = reverse(ind) + 1;
                    %disp('reverse')
                    %disp(N(j))
            end
            %m_dot(i) = 0;
            
        end   

        d_mass_d_time = m_dot(i);
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
    veff2(ind) =  mass(end)/mass(1);
    veff1(ind) = (1 - mass_start/mass(1));
    V_eff(ind) = V_eff(ind) * veff1(ind) * veff2(ind);
    
    ind = ind + 1;
    %m_dot * temp * dt + temp in * mass
end
