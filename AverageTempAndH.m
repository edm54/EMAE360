%Returns the cumulative work for one cylinder in kJ, given a (volumetric
%efficiency, equivalence ratio, RPM).
function qavg =  AverageTempAndH(nv,eq,N,plt)
% Real Engine Cycle with Heat Transfer from the Cylinder
step = 1; %crankangle interval for calculation
r=10;
thetai = 0; %initial crankangle, rad
thetaf =721; %final crankangle througout steps
NN = (thetaf - thetai)/step + 1;

thetas = 160; %start of heat release (deg)
thetad = 60; %duration of heat release (deg)

af = 14.7/eq;

ma = 0.000336*nv; %kg
mf = ma/af;
b = 0.0705; %Bore m
s = 0.0641; %Stroke m
Ao = 1.5*pi*(b/2)^2; %Total Area of Combustion chamber at TDC, will come from 
Vd = 0.0015/6;
V1 = r/(r-1)*Vd;
Vo = V1/r; %Total volume at top dead center, m^3
P1 = 101.325; %Initial Pressure, kPa
T1 = 300; %K
Tw = 530; %K, inside wall temperature
Tbar = Tw / T1;
Qlhv = 44000000; %J/kg
nc = 0.95;
Qin = mf*Qlhv*nc/(1000*P1*V1); %Dimensionless total Heat release
omega = N*pi/60; %rps
Ubar = 2*s*N/60; %Mean piston speed (m/s)
B = 4*r/(b*(Ao/Vo)-4); %dimensionless volume

prop.theta = thetai:step:thetaf;
prop.theta = prop.theta.';
prop.vol = zeros(NN,1);
prop.area = zeros(NN,1);
prop.press = zeros(NN,1);
prop.temp = zeros(NN,1);
prop.rho = zeros(NN,1);
prop.work = zeros(NN,1);
prop.heatloss = zeros(NN,1);
prop.htcoeff = zeros(NN,1);
prop.heatflux = zeros(NN,1);
gam = zeros(NN,1);

prop.vol(1) = nondimV(0,r);
prop.area(1) = (Ao-4*Vo/b) + 4* prop.vol(1)/b;
prop.press(1) = 1; %nondimensional initial pressure, P1/P1
prop.temp(1) = T1; %initial temperature, K

fy = [prop.press(1);0;0]; %Initial Pressure, Initial work, heat loss

for i=1:721
    [prop.vol(i),~]=nondimV(prop.theta(i),r); 
end

for i = 1:360
    [fy,ht,hflux,prop.temp(i)] = integrate_ht(prop.theta(i),prop.theta(i)+step,fy);
    gam(i) = calc_gamma(prop.temp(i));
    prop.press(i)=fy(1);
    prop.work(i)=fy(2);
    prop.heatloss(i)=fy(3);
    prop.htcoeff(i)=ht;
    prop.heatflux(i)=hflux;
end

U = 6.18*Ubar;
for i=360:1:540
   prop.temp(i)=prop.temp(360);
   prop.htcoeff(i) = heatTransferCoeff(prop.press(360)*P1,U,b,prop.temp(i));
   prop.heatflux(i) = prop.htcoeff(i)*T1*(prop.press(360)*prop.vol(i) - Tbar)/10^6;
end

for i=540:1:721
   prop.temp(i)=300;
   prop.htcoeff(i)=heatTransferCoeff(101,U,b,300);
   prop.heatflux(i) = prop.htcoeff(i)*T1*(1*prop.vol(i) - Tbar)/10^6;
end
h_avg = zeros(NN,1);
T_avg = zeros(NN,1);
q_avg = zeros(NN,1);

prod = prop.temp.*prop.htcoeff;
havg = (1/(720))*trapz(prop.htcoeff);
h_avg(:,1) = havg;
Tavg = (1/(720*havg))*trapz(prod);
T_avg(:,1) = Tavg;
qavg = (1/(720))*trapz(prop.heatflux);
q_avg(:,1) = qavg;

if plt>0
    figure
    plot(prop.theta,prop.temp,prop.theta,T_avg)
    ylabel("Temperature, K")
    xlabel("Angle, degree")
    title("Temperature vs Crank Angle")
    legend("Instantaneous Temperature", "Average Temperature")


    figure
    plot(prop.theta,prop.heatflux,prop.theta,q_avg)
    ylabel("q'', MW/m^2")
    xlabel("Angle, degree")
    title("Heat Flux vs Crank Angle")
    legend("Instantaneous Heat Flux", "Average Heat Flux")

    figure
    plot(prop.theta,prop.htcoeff,prop.theta,h_avg)
    ylabel("hg, W/(m^2 K)")
    xlabel("Angle, degree")
    title("Gas side Heat transfer coefficient vs Crank Angle")
    legend("Instananeous hg", "Average hg")
    
    press = prop.press;
    press = press .* 0.001 .* P1;
    figure
    plot(prop.theta, press)
    xlabel("Angle, deg")
    ylabel("Pressure, MPa")
    title("Pressure vs Crank Angle")
end

    function [fy,ht,hflux,T] = integrate_ht(theta,thetae,fy)
        [ht,hflux,T] = values(theta,fy);
        [tt,yy] = ode23(@rates,[theta thetae],fy);
        for j=1:3
           fy(j) = yy(length(tt),j); 
        end
    end

    function [yprime] = rates(theta, fy)

        [V,dV_dtheta] = nondimV(theta,r);
        P = P1*fy(1);
        T = T1*fy(1)*V;
        gamma = calc_gamma(T);
        U = 2.28 * Ubar + 0.00324*T1*(r-1)*(fy(1)-V^(-gamma))/r;
        dxb_dtheta = 0;
        if theta > thetas 
            [~, dxb_dtheta] = burnRate(theta, thetas, thetad);
        end
        ht = heatTransferCoeff(P,U,b,T);
        h = ht*T1*4/(P1*1000*omega*B*b);

        term1= -gamma*fy(1)*dV_dtheta/V;
        term3= h*(1. + B*V)*(fy(1)*V - Tbar)*pi/180.;
        term2= (gamma-1)/V*(Qin*dxb_dtheta - term3);
        yprime(1,1)= term1 + term2 ;
        yprime(2,1)= fy(1)*dV_dtheta;
        yprime(3,1)= term3;
    end

    function [ht,hflux,T] = values(theta, fy)
        [V,~] = nondimV(theta,r);
        P = P1*fy(1);
        T = T1*fy(1)*V;
        gamma = calc_gamma(T);
        U = 2.28 * Ubar + 0.00324*T1*(r-1)*(fy(1)-V^(-gamma))/r;
        ht = heatTransferCoeff(P,U,b,T);
        hflux=ht*T1*(fy(1)*V - Tbar)/10^6;
    end

    function h = heatTransferCoeff(P,U,b,T)
    h = 3.26*(P^0.8)*(U^0.8)*(b^(-.2))*(T^(-.55));
    end

    function [V,rate] = nondimV(theta,r)
        V = (1/r) + (r-1)/(2*r).*(1+cosd(theta));
        rate = -(r-1)/(2*r)*sind(theta)*pi/180;
    end

    function [xb, rate] = burnRate(theta, thetas, thetad)
    a = 5;
    n = 3;
    xb = 1-exp((-a*((theta - thetas)/thetad)^n));
    rate = n*a*(1-xb)/thetad*((theta - thetas)/thetad)^(n-1);
    end
end