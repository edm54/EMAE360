function [volume] =  calc_volume(theta)
bore = 7.05;
r_bs = 1.1;
stroke = bore/r_bs;
crank_rad = stroke/2; %cm
% C rad/length ~ = 1/3
crank_l = 3.5 * crank_rad; %cm

rc = 10;

% Vc = Vtotal / (rc -1)
Vc = (1500/6)/(rc -1);

s = crank_rad * cos(theta*pi/180) + sqrt(crank_l^2 + (crank_rad^2) * (sin(theta*pi/180))^2); 
volume = (Vc + (crank_l + crank_rad - s)*(pi*bore^2)/4)/1e6;

end
