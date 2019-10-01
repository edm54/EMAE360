%%
% This program calculates bore, stroke based on a given ratio
% Also calculates crank radius and length 
% Calculates volume as a function of time

clear all
clc
close all

total_displacement = 1500
%for total_d = 1500: 25: 1800
S = 1500 %cM/sec
N = 5000/60; %rev/sec

% Heywood 44
stroke = S/(2 * N) %cM

displacement = total_displacement / 6 %cubic centimeters

% Heywood p 51, d = b^2* pi * L/4
bore = sqrt(4 * displacement /(pi*stroke)) %cm
R_bs = bore/stroke %Unitless bore stroke ratio 

% orrr
bs_ratio = 1.1
stroke = (sqrt(4*displacement)/(sqrt(pi)*bs_ratio))^(2/3)
mean_p_speed = stroke * 2 * 5000/60
max_p_speed = stroke * 2 * 9400/60
bore = sqrt(4 * displacement /(pi*stroke))

%%
% Heywood p 43
crank_rad = stroke/2 %cm
% C rad/length ~ = 1/3
crank_l = 3.5 * crank_rad %cm

rc = 10

% Vc = Vtotal / (rc -1)
Vc = (1500/6)/(rc -1)

i = 1
for theta = -180:1:180
    % Heywood 44
    s(i) = crank_rad * cos(theta*pi/180) + sqrt(crank_l^2 + (crank_rad^2) * (sin(theta*pi/180))^2); 
    % Heywood 43
    v(i) = Vc + (crank_l + crank_rad - s(i))*(pi*bore(1)^2)/4;
    i = i + 1;
end


%%

%45 of heywood: 
N = 9000
mean_piston_speed = 2 * stroke/100 * N/60
R = 3.5 
i = 1
for theta = -180:1:180
    % Heywood 44
    piston_speed(i) = mean_piston_speed * (pi/2) * sin(theta*pi/180)*(1+ ((cos(theta*pi/180))/sqrt(R^2 - (sin(theta*pi/180))^2)))
    i = i + 1;
end
figure
% TEmp asa fun of angle
plot(0:1:360, piston_speed)
xlim([0,360])
xlabel('Crank Angle (Degrees)')
ylabel('Instantaneous Piston Velocity (m/s)')
title('Instantaneous Piston Velocity vs Crank Angle')

figure

%%
plot(R_bs, total_displacement : 25 : 1800)
xlabel('Bore-Stroke Ratio')
ylabel('Displacement (cc)')
title('Bore Stroke Ratio as a Function of Displacement')
figure
 
% TEmp asa fun of angle
plot(0:1:360, v)
xlim([0,360])
xlabel('Crank Angle (Degrees)')
ylabel('Volume of Combustion Chamber (cc)')
title('Volume of Cylinder vs Crank Angle')
st = strcat('Minimum volume =  ' , num2str(v(1)), ' cc')
ht = text(110, v(1), st);
st = strcat('Maximum volume =  ' , num2str(v(181)), ' cc')
ht = text(110, v(181)+10, st);