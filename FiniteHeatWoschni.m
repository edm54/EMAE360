%Returns the cumulative work for one cylinder in kJ, given a (volumetric
%efficiency, equivalence ratio, RPM).
function Wnet = FiniteHeatWoschni(Qin,N,ma,plt)

%Motored Engine
r = 10;
T1 = 300; %K 
P1 = 101325; %Pa

% Real Engine Cycle with Heat Transfer from the Cylinder
step = 1; %crankangle interval for calculation
thetai = 0; %initial crankangle, rad
thetaf =360; %final crankangle througout steps
NN = (thetaf - thetai)/step + 1;

thetas = 160; %start of heat release (deg)
thetad = 60; %duration of heat release (deg)

b = 0.0705; %Bore m
s = 0.0641; %Stroke m
Ao = 3850/(1000)^2; %Total Area of Combustion chamber at TDC, m
Vd = 0.0015/6;
V1 = r/(r-1)*Vd;
Vo = V1/r; %Total volume at top dead center, m^3
P1 = 101.325; %Initial Pressure, kPa
T1 = 300; %K
Tw = 530; %K, inside wall temperature
Tbar = Tw / T1;
omega = N.*pi/60; %rps
Ubar = 2.*s.*N./60; %Mean piston speed (m/s)
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

for i = 2:NN
    [prop.vol(i),~] = nondimV(prop.theta(i), r);
    prop.rho(i) = ma / (V1*prop.vol(i));
    prop.area(i) = (Ao-4*Vo/b) + 4* prop.vol(i)*V1/b;
end

for i = 1:NN
    [fy,ht,hflux,prop.temp(i)] = integrate_ht(prop.theta(i),prop.theta(i)+step,fy);
    gam(i) = calc_gamma(prop.temp(i));
    prop.press(i)=fy(1);
    prop.work(i)=fy(2);
    prop.heatloss(i)=fy(3);
    prop.htcoeff(i)=ht;
    prop.heatflux(i)=hflux;
end

bar = prop.press .* P1/100; %Bar
prop.press = prop.press.*P1/(10^3); %MPa
work = P1.*V1 .* prop.work; %kJ
heatloss = P1.*V1.*prop.heatloss; %kJ
Wnet = work(NN);

if plt>0
    thetai = 0;
    mstep = 0.001;
    thetaf = 360;
    mprops.theta = thetai:mstep:thetaf;
    mprops.theta = mprops.theta.';
    mprops.vol = zeros(length(mprops.theta),1);
    mprops.press = zeros(length(mprops.theta),1);
    mprops.temp = zeros(length(mprops.theta),1);
    mprops.rho = zeros(length(mprops.theta),1);
    mgamma = zeros(length(mprops.theta),1);

    mprops.vol(1) = nondimV(thetai, r)*V1;
    mprops.press(1) = 1;
    mprops.temp(1) = T1;
    k = P1*(V1)^(calc_gamma(T1));
    mprops.rho(1) = P1/(287*T1);
    for i = 2:length(mprops.theta)
        mgamma(i) = calc_gamma(mprops.temp(i-1));
        [mprops.vol(i),~] = nondimV(mprops.theta(i),r);
        mprops.vol(i) = mprops.vol(i)*V1;
        mprops.rho(i) = ma / (mprops.vol(i));
        mprops.press(i) = k / (mprops.vol(i))^(mgamma(i));
        mprops.temp(i) = mprops.press(i)/(mprops.rho(i)*287);
    end
    
    mprops.press = mprops.press ./ (10^6);
    figure
    hold on
    plot(prop.theta,prop.press)
    plot(mprops.theta,mprops.press)
    ylabel("Pressure, MPa")
    xlabel("Angle, rad")
    title("Pressure vs Crank Angle")
    legend("Actual Pressure", "Motored Pressure")
    hold off 

    figure
    plot(prop.theta,bar)
    ylabel("Pressure, bar")
    xlabel("Angle, deg")
    title("Pressure vs Crank Angle")

    figure
    plot(prop.theta,prop.temp)
    ylabel("Temperature, K")
    xlabel("Angle, deg")
    title("Temperature vs Crank Angle")
    %{
    figure
    plot(prop.theta,gam)
    ylabel("Gamma")
    xlabel("Angle, rad")
    title("Gamma vs Crank Angle")
    %}

    figure
    plot(prop.theta,prop.heatflux)
    ylabel("q'', MW/m^2")
    xlabel("Angle, deg")
    title("Heat loss vs Crank Angle")

    figure
    plot(prop.theta,prop.htcoeff)
    ylabel("h")
    xlabel("Angle, deg")
    title("Heat transfer coefficient vs Crank Angle")

    figure
    plot(prop.theta,work,prop.theta,heatloss)
    legend('Work','Heatloss')
    xlabel('Crank Angle, deg')
    ylabel('Cumulative Work and Heat Loss per Cylinder, kJ')
    title('Cumulative Work and Heat Loss vs Crank Angle per Cylinder')
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