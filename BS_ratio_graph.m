clear all
set(0, 'DefaultAxesFontWeight', 'normal', ...
      'DefaultAxesFontSize', 18, ...
      'DefaultAxesFontAngle', 'normal', ... 
      'DefaultAxesFontWeight', 'normal', ... 
      'DefaultAxesTitleFontWeight', 'bold', ...
      'DefaultAxesTitleFontSizeMultiplier', 1.2) ;
set(groot,'defaultLineLineWidth',3)

displacement = 1500/6
%bs_ratio = 1.1;
figure
for bs_ratio = .7 : .2 : 1.5
    i = 1;
    for N = 2000 : 100 : 12000
        stroke = (sqrt(4*displacement)/(sqrt(pi)*bs_ratio))^(2/3);
        S(i) = (stroke * 2 * N/60)/100;
        i = i + 1 ;
    end
    if bs_ratio == 1.1
        plot(2000 : 100 : 12000, S, 'LineWidth', 5)
    else
        plot(2000 : 100 : 12000, S)
    end
    
    hold on
    
end
ylim([0 25])
xlabel('Revolutions per minute (RPM)')
ylabel('Mean Piston Speed (M/S)')
legend('.7', '.9', '1.1', '1.3', '1.5')
title('Mean Piston Speed vs RPM for Varying Bore/Stroke Ratio')
x=2000:10:12000;
y=20;
plot(x,y*ones(size(x)), '--k','LineWidth', 2, 'HandleVisibility','off')

%line([2000,12000],[20, 20] , 'LineWidth', 2, 'HandleVisibility','off')


