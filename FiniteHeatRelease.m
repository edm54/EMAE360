function[v]=FiniteHeatRelease()
% Gas cycle heat release code for two engines
close all
% engine parameters 
clear();
%%
calc_gamma(300)
%%


thetas(1,1)= -10; % Engine1 start of heat release (deg)
thetas(2,1)= -10;   % Engine2 start of heat release (deg)
thetad(1,1) = 40; % Engine1 duration of heat release (deg)
thetad(2,1) = 10; % Engine2 duration of heat release (deg)

r=10;       %compression ratio
gamma= 1.4; %gas const
q= 34.385233793707900;      % dimensionless total heat release Qin/P1V1 

a= 5;  %weibe parameter a 
n= 3; %weibe exponent n

step=1;      % crankangle interval for calculation/plot
NN=360/step; % number of data points

% initialize the results data structure
save.theta=zeros(NN,1); % crankangle
save.vol=zeros(NN,1);   % volume 
save.press=zeros(NN,2); % pressure 
save.work=zeros(NN,2);  % work
save.temp = zeros(NN,2);


pinit(1) = 1; % Engine 1 initial dimensionless pressure P/P1
pinit(2) = 1; % Engine 2 initial dimensionless pressure P/P1 
ma = 3.335683867042752e-04;
n_constant = ma/(1000*28.97);
r_constant = 8.314;
stroke = 6.4100000000000; %cm
bore = 7.051000000000001; %cm
bore = 1.1 * stroke;
crank_rad = stroke/2; %cm
crank_l = 3.5 * crank_rad; %cm
Vc = (1500/6)/9;

i = 1;
for theta = -180:1:180
    s(i) = crank_rad * cos(theta*pi/180) + sqrt(crank_l^2 + (crank_rad^2) * (sin(theta*pi/180))^2); 
    v(i) = (Vc + (crank_l + crank_rad - s(i))*(pi*bore(1)^2)/4)/1e6;
    i = i + 1;
end

% for loop for engine1 and engine2
for j=1:2
theta = -180;          %initial crankangle
thetae = theta + step; %final crankangle in step
fy(1) = pinit(j); % assign initial pressure to working vector
fy(2) =0.;        % reset work vector 

% for loop for pressure and work calculation
  for i=1:NN
    [fy, vol] = integrate(theta,thetae,fy);
    
    % reset to next interval
    theta = thetae;
    thetae = theta+step;
    % copy results to output vectors
    save.theta(i)=theta;
    save.vol(i)=vol;
    save.press(i,j)=fy(1);
    save.work(i,j)=fy(2);
    %save.temp(i,j) = save.press(i,j) * v(i)/(n_constant * r_constant);
    P = 1e5 * save.press(i,j);
    R = 287.058;
    rho(i) = ma/(v(i));
    save.temp(i,j) = P/(R*rho(i));
    disp(rho(i))
    gamma = calc_gamma(save.temp(i,j));
    g(i) = gamma;
    
  end %end of pressure and work iteration loop
end %end of engine iteration loop

[pmax1, id_max1] = max(save.press(:,1)); %Engine 1 max pressure
[pmax2, id_max2] = max(save.press(:,2)); %Engine 2 max pressure
thmax1=save.theta(id_max1);%Engine 1 crank angle
thmax2=save.theta(id_max2);%Engine 2 crank angle

w1=save.work(NN,1);
w2=save.work(NN,2);
eta1= w1/q; % thermal efficiency
eta2= w2/q; 
imep1 = eta1*q*(r/(r -1)); %imep
imep2 = eta2*q*(r/(r -1));
eta_rat1 = eta1/(1-r^(1-gamma)); 
eta_rat2 = eta2/(1-r^(1-gamma));



% output overall results
fprintf('                  Engine 1           \n');
fprintf(' Theta_start       %5.2f              \n', thetas(1,1));
fprintf(' Theta_dur          %5.2f                \n', thetad(1,1));
fprintf(' P_max/P_1           %5.2f              \n', pmax1);
fprintf(' Theta_max       %7.1f             \n',thmax1);
fprintf(' Net Work/P1V1   %7.2f             \n', w1);
fprintf(' Efficiency         %5.3f               \n', eta1);
fprintf(' Eff. Ratio         %5.3f              \n', eta_rat1 );
fprintf(' Imep/P1            %5.2f                \n', imep1);

%plot results

%set(gcf,'Units','pixels','Position', [50,50,1200,600]);
%subplot(1,2,1);
plot(save.theta,save.press(:,1),'linewidth',3 )
set(gca, 'fontsize', 18,'linewidth',2);
%grid
legend('Engine 1','Location','NorthWest')
xlabel('Theta (deg)','fontsize', 18)
ylabel('Pressure (bar)','fontsize', 18)
title('Pressure (with changing gamma) as a Function of Crank Angle')
grid
print -deps2 heatrelpressure;

figure();
%subplot(1,2,2);
plot(save.theta,save.work(:,1),'-', 'linewidth',2)
set(gca, 'fontsize', 18,'linewidth',2);
%grid
legend('Engine 1', 'Location','NorthWest')
xlabel('Theta (deg)','fontsize', 18)
ylabel('Work','fontsize', 18)
title('Work as a Function of Crank Angle')
grid

figure();
%subplot(1,2,2);
plot(save.theta,save.temp(:,1),'-', 'linewidth',2)
set(gca, 'fontsize', 18,'linewidth',2);
%grid
legend('Engine 1', 'Location','NorthWest')
xlabel('Theta (deg)','fontsize', 18)
ylabel('Temp','fontsize', 18)
title('Temperature (with changing gamma) as a Function of Crank Angle')
grid

figure();
%subplot(1,2,2);
plot(save.theta,save.temp(:,1),'-', 'linewidth',2)
set(gca, 'fontsize', 18,'linewidth',2);
%grid
legend('Engine 1', 'Location','NorthWest')
xlabel('Theta (deg)','fontsize', 18)
ylabel('Temp','fontsize', 18)
title('Temperature (with changing gamma) as a Function of Crank Angle')
grid

figure();
%subplot(1,2,2);
plot(save.theta,g,'-', 'linewidth',3)
set(gca, 'fontsize', 18,'linewidth',2);
%grid
%legend('Engine 1', 'Location','NorthWest')
xlabel('Theta (deg)','fontsize', 18)
ylabel('Temp','fontsize', 18)
title('Gamma as a Function of Crank Angle')
grid


function[fy,vol] = integrate(theta,thetae,fy)
%  ode23 integration of the pressure differential equation 
%  from theta to thetae with current values of fy as initial conditions
[tt, yy] = ode23(@rates, [theta thetae], fy);
%put last element of yy into fy vector
 for k=1:2
  fy(k) = yy(length(tt),k);
 end
%pressure differential equation
    function [yprime] = rates(theta,fy) 
    vol=(1.+ (r -1)/2.*(1-cosd(theta)))/r;
    dvol=(r - 1)/2.*sind(theta)/r*pi/180.; %dvol/dtheta
    dx=0.; %set heat release to zero
        if(theta > thetas(j)) % then heat release dx  > 0
        dum1=(theta -thetas(j))/thetad(j);
        x=1.- exp(-(a*dum1^n));
        dx=(1-x)*a*n*dum1^(n-1)/thetad(j); %dx/dthetha
        end  
     term1= -gamma*fy(1)*dvol/vol;
     term2= (gamma-1)*q*dx/vol;
     yprime(1,1)= term1 + term2;
     yprime(2,1)= fy(1)*dvol;    
    end %end of function rates
end  %end of function integrate2


function [gamma] =  calc_gamma(temp)
    N2_high = [.2896e1, .15155e-2, -.57235e-6, .99807e-10, -.652e-14];
    N2_low = [.3675e1, -.12082e-2, .23240e-5, -.6322e-9, -.2258e-12];
    O2_high = [.362e1, .7362e-3, -.1965e-6, .362e-10, -.2895e-14];
    O2_low = [.36256e1, -.18782e-2, .70555e-5, -.6764e-8, .21556e-11];

    if (temp>1000)
        N = N2_high;
        O = O2_high;
    else
        N = N2_low;
        O = O2_low;
    end    
    
    cp_n = temp_to_gamma(N, temp);
    cp_o = temp_to_gamma(O,temp);
    
    cp = .2095 * cp_o + .7905 * cp_n;
    R = .287;
    cp = cp*R;
    cv = cp - R;
    function [cp] = temp_to_gamma(a, temp)
        cp =  a(1) + a(2)*temp + a(3)*temp^2 + a(4)*temp^3 + a(5) * temp^4;
    end
    gamma = cp/cv;

end
end % heat_release_weibe2

