%Air intate
clear all
clc
close all

%% Prepare Pllots
set(0, 'DefaultAxesFontWeight', 'normal', ...
      'DefaultAxesFontSize', 18, ...
      'DefaultAxesFontAngle', 'normal', ... 
      'DefaultAxesFontWeight', 'normal', ... 
      'DefaultAxesTitleFontWeight', 'bold', ...
      'DefaultAxesTitleFontSizeMultiplier', 1.2) ;
set(groot,'defaultLineLineWidth',3)


figure
hold on
l_in1 = [.05 .08 .1];
c_1 = [.58 .62 .67];

l_in2 = [.1 .11];
c_2 = [.67 .55] ;

l_in3 = [.11 .14 .18];
c_3 = [.55 .59 .68];

l_in4 = [.18 .21 .24 .28];
c_4 = [.68 .55 .53 .48];

fit1 = polyfit(l_in1, c_1,  3);
val1 = polyval(fit1, .05:.01:.1);
% plot(.05:.01:.1, val1);
% hold on
fit2 = polyfit(l_in2, c_2,  1);
val2 = polyval(fit2, l_in2);
% plot(l_in2, val2);

fit3 = polyfit(l_in3, c_3,  2);
val3 = polyval(fit3, .11:.01:.18);
% plot(.11:.01:.18, val3);

fit4 = polyfit(l_in4, c_4,  1);
val4 = polyval(fit4, .18:.01:.35);
% plot(.18:.01:.35, val4);


bore = 7.045/100;
r_bs =  1.1;
stroke = bore/r_bs;


opt1 =  .00719;
opt2 = .00733;

ex_mult = linspace(.005, .0085, 5)
ex_mult = [.006 .007  .008]
for index1 = 1:length(ex_mult)
    %max_in_lift = ex_mult(index1) * bore; %meter
    max_in_lift = .007 %meter
    
    max_in_lift = ex_mult(index1)
    max_ex_lift = .008459 % meter
  
    %%max_v_lift = .10 * bore; 
    %max_v_lift = .08 * bore;

    theta = [-360 : .25: 360];
    y1 = -1*(theta.^2 -220*theta - 1125);

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

%     plot(theta, exhuast_lift*1000);
%     plot(theta,intake_lift*1000);
%     title('Intake Valve Lift as a function of crank angle')
%     ylim([0,10]);
%     legend('Exhuast valve lift','Intake valve lift')
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    %saveas(gcf,'Valve Lift vs Crank.png')

    %%
    % Look at 5000 hieght and AOC 
    % Max areas:
       % In: 2.6 cm
       % Ex: 2.25 cm
       % Min is 1 cm 
       % .8459 cm of lift
        % floor of .1 cm
    ex_mult2 = [.023 .026]
    for index2 = 1:length(ex_mult2)
        D_in = .026;
        D_in = ex_mult2(index2);
        D_ex = .0225 ;

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


    %     plot(theta,Am_in)
    %     plot(theta, Am_ex)
    % 
    %     legend('Intake', 'Exhuast')
    %     title('Flow Area (m^2) as a Function of Crank Angle')
    %     legend('Intake Flow Area', 'Exhuast Flow Area')
    %     ylabel('Flow Area (M^2)')
    %     xlabel('Crank Angle')
    %     %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    %     %saveas(gcf,'Flow Area vs Crank.png')
    %     figure

        %% ex
        cd_ex = .5 + .5 .* [  1.2 1.5 1.6 1 .4 .2 0]./3.2;


        lv_dv_ex = [.1 .15 .2 .25 .3 .35 .4];
    %     plot(lv_dv_ex, cd_ex);
        fit = polyfit(lv_dv_ex, cd_ex,  3);
        val = polyval(fit, .1:.01:.3);
    %     hold on
    %     plot(.1:.01:.3, val);
    %     ylim([0 1]);
    %     xlim([0 .3]);

        %% intake Cd vs lv/dv
        %figure
        %title('Intake Cd vs LV/DV')
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
    %     figure
    %     plot(theta, CD_in)
    %     hold on
    %     plot(theta, CD_ex)
    %     title('CD as a Function of Crank Angle')
    %     xlabel('Crank angle')
    %     ylabel("Discharge Coeff")
    %     legend('Intake', 'Exhuast')
        %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
        %saveas(gcf,'CD vs Crank.png')

        %% Exhuast
        ind = 1;
        %N = 800: 25: 10000;
        %N = 800:25:10000;
        N = 800:100:9400;
        %N = 5000
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

            %figure
            %hold on
            plot(N, V_eff)
    %         plot(N, veff1)
    %         plot(N, veff2)
    %         grid 
    %         s = sprintf('Volumetric Efficency vs RPM for lift = %.5f', max_ex_lift);  
    %         title(s)
    %         legend('Total', 'Exhuast', 'Intake')
            % This legend is for Exhuast lift
            %legend('0.00423','0.0053','0.0063','0.0074','0.00845')
            [max_v_eff(index1), max_ind(index1)] = max(V_eff)
            max_v_eff_ex(index1) = max(veff1)
            max_v_eff_in(index1) = mean(veff2)
            m_lift_ex(index1) = max_in_lift
            AOC(index1) = trapz(N,V_eff)
    %        disp(ex_mult(index1))
    %        disp(max(V_eff))
    %        disp(max(veff1))
    %        disp(max(veff2))
         
           legend('.005',  '.006', '.007', '.008', '.009')
           legend('.0045', '.0055', '.0065', '.0075',   '.0085')
           legend('.023 .006', '.026, .006', ...
                  '.023, .007', '.026, .007',...
                  '.023, .008', '.026, .008')
           title('Volumetric Efficiency vs Exhuast Valve Lift')
           xlabel('RPM')
           ylabel('Volumetric Efficiency')
           %s = sprintf('For max_lift = %.5f, V total = %.5f, V Ex = %.5f, V In = %.5f', max_ex_lift, max(V_eff), max(veff1), max(veff2))  
           ex_mult2(index2)
           ex_mult(index1)
    end
    
end
%%
figure
%plot(m_lift_ex, max_v_eff_in)
hold on
%plot(m_lift_ex, max_v_eff_ex)
plot(m_lift_ex, AOC)
%plot(m_lift_ex, N(max_ind))
%plot(m_lift_ex, max_v_eff)
%legend('In', 'ex', 'total')
title('Area Under Volumetric Efficiency Curve vs Exhuast Lift')


%%
figure 
plot(theta,Am_in .* CD_in) 
hold on
plot(theta,Am_ex .* CD_ex) 
title('Min Area * Discharge Coefficient vs. Crank Angle')
