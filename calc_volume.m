function [volume] =  calc_volume(theta)

bore = .0705
s = crank_rad * cos(theta*pi/180) + sqrt(crank_l^2 + (crank_rad^2) * (sin(theta*pi/180))^2); 
volume = (Vc + (crank_l + crank_rad - s)*(pi*bore^2)/4)/1e6;

end
