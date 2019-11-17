%{
This function returns the net work required to operate a "turned off"
cylinder for 1 revolution. 
%}
function Wnet = PumpingLoss(ma)
    ho = 426.4; %kJ/kg
    Wnet = ma*ho;
end