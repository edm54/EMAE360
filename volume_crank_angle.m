function [volume] = volume_crank_angle(theta)
%Returns the volume of the chamber given a crank angle in degrees

Vc = (1500/6)/(rc -1)

% Heywood 44
s = crank_rad * cos(theta*pi/180) + sqrt(crank_l^2 + (crank_rad^2) * (sin(theta*pi/180))^2); 
% Heywood 43
volume = Vc + (crank_l + crank_rad - s(i))*(pi*bore(1)^2)/4;


end

