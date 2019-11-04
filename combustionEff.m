function [combustion_eff] = combustionEff(eq_ratio)
%UNTITLED4 Summary of this function goes here
%   Recreation of Figure 4.1 from Pulkabrek 
% Fuel Efficicency ratio plot

x = [1 1.1 1.2 1.3 1.4]
y = [ .97 .9 .8 .7 .6]
fit = polyfit(x, y,  2);
val = polyval(fit, 1:.01:1.4)

if eq_ratio>1
    combustion_eff = polyval(fit, eq_ratio)
else
    combustion_eff = .97
end
    
end