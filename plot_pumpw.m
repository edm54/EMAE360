
N = 800:100:9400;
P_atm = 101325;
%[veff, pump] = volumetric_efficiency(N)
pressure_drops = [101000 70000  40000 15000]



for i = 1:length(pressure_drops)
    [veff, pump] = volumetric_efficiency_p(N, pressure_drops(i))
%     for j = 1:length(N) 
%         vol(i,j) = veff(j);
%         pw(i,j) = pump(j);
%         php(i,j) = (pump(j)/1000 * N(i)/120)/0.745699872;   
     php(i,:) = 6  * (pump.'/1000 .* N/120)/0.745699872   
end

equivolence_ratio = 1
combustion_eff = combustionEff(equivolence_ratio)

RPM = [2100 15000];
mech_eff = [ .9 .65];  
for i = 1:length(N)  
    c(i) = combustion_eff
    if N(i)<= 2100
        mechanical_eff(i) = .9;
    else
        mechanical_eff(i) = interp1(RPM, mech_eff, N(i), 'linear');
    end
end
%%
figure
plot(N, mechanical_eff)
grid
hold on 
plot(N, veff)
plot(N,c)
plot(N, min(1, pl))
title('Engine Efficiency Losses')
xlabel('RPM')
ylabel('Efficiency')
legend('Mechanical', 'Volumetric','Combustion', 'Pumping Loss') 

%%
for iN = 1:length(pressure_drops)
      legendCell{iN} = num2str(pressure_drops(iN),'P intake=%-d');
end

figure
plot(N, php)
legend(legendCell)
title('Pumping Loss (hp) Varying Intake Pressure')
xlabel('RPM')
ylabel('Power Lost to Pumping (hp)')


