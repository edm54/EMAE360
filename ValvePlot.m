clear all
%% EXHUAST
theta = [-360 : .1: 360];
b = 0
h = 8.5
t0 = 0
t1 = pi/4
intake = zeros(1, length(theta))
exhuast = zeros(1, length(theta))
% t0 = -220

t0 = -260
t1 = -105
tend = 50
half = (t0 + t1)/2

B = t1 - t0
i = 1
exhuast = zeros(1,length(theta))
%%
%for x = 0:.01:pi/4 + .01
for x = t0:.1:t1
    R(i) = b + (h/(2 * pi)) * ((2 * pi/ B)*( x - t0) - sin((2*pi / B) * ( x - t0)));
    i = i+1
end
j = i 
figure
%plot(-220:1:-105, R)
plot(t0:.1:t1, R)

i = 1
for x = t0:.1:t1
    R1(i) =  h +  b - (h/(2 * pi)) * ((2 * pi/ B)*( x - t0) - sin((2*pi / B) * ( x - t0)));
    i = i+1
end

i = 1
for j = 1:length(theta)
    if theta(j)>= t0 && theta(j) <= t1
        exhuast(j) = R(i);
        i = i + 1;
    end
    
    if theta(j) == t1
        i = 1;
    end
    
    if theta(j)> t1 && theta(j)<tend
        exhuast(j) = R1(i);
        i = i +1;
    end
    
end


hold on
plot(-105:.1:50, R1)

%% Intake

b = 0
h = 7
t0 = 0
t1 = pi/4

t0 = -5
t1 = 215
tend = 245
t0 = -35
t1 = 105
half = (t0 + t1)/2

B = t1 - t0
i = 1
%for x = 0:.01:pi/4 + .01
for x = t0:.1:t1
    R2(i) = b + (h/(2 * pi)) * ((2 * pi/ B)*( x - t0) - sin((2*pi / B) * ( x - t0)));
    i = i+1
end

%figure
plot(t0:.1:t1, R2)
j = i
i = 1
for x = t0:.1:t1
    R3(i) =  h +  b - (h/(2 * pi)) * ((2 * pi/ B)*( x - t0) - sin((2*pi / B) * ( x - t0)));
    i = i+1
end

intake = zeros(1,length(theta))
i = 1
for j = 1:length(theta)
    if theta(j)>= t0 && theta(j) <= t1
        intake(j) = R2(i);
        i = i + 1;
    end
    
    if theta(j) == t1
        i = 1;
    end
    
    if theta(j)> t1 && theta(j)<tend
        intake(j) = R3(i);
        i = i +1;
    end
end

hold on
plot(t1:.1:245, R3)

title('Smoothed Lift Curve')
xlabel('Lift (mm)')
ylabel('Crank Angle')
legend('Exhuast', 'Exhuast', 'Intake', 'Intake')
%%
bore = 7.045/100;
r_bs =  1.1;
stroke = bore/r_bs;
%max_in_lift = .12 * bore; %meter
%max_ex_lift = .12 * bore;
max_in_lift = .007 %meter
max_ex_lift = .008459 % meter

theta = [-360 : .1: 360];
y1 = -1*(theta.^2 -220*theta - 1125);

open_in = -5
in_dur = 220

close_in = in_dur + open_in
intake_lift = -1 * (theta - open_in).* (theta - close_in); 

open_ex = -220
ex_dur = 230
close_ex = ex_dur + open_ex;

exhuast_lift = -1 * (theta - open_ex).* (theta - close_ex); 
y_out_max = max(exhuast_lift)
y_in_max = max(intake_lift)

intake_lift = (max_in_lift/y_in_max).* intake_lift
exhuast_lift = (max_ex_lift/y_out_max) .* exhuast_lift

for i = 1:length(intake_lift)
    if intake_lift(i) < 0
        intake_lift(i) = 0;
    end
    if exhuast_lift(i) < 0
        exhuast_lift(i) = 0;
    end
end

%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%saveas(gcf,'Valve Lift vs Crank.png')
figure
plot(theta, exhuast_lift*1000);
hold on
plot(theta,intake_lift*1000);

%%
figure 
plot(theta, exhuast_lift*1000)
hold on
% plot(-105:.1:50, R1)
% plot(-260:.1:-105, R)
plot(theta, exhuast)
%%
figure 
plot(theta, intake_lift*1000)
hold on
plot(t0:.1:t1, R2)
plot(t1:.1:245, R3)

%%
figure
plot(theta, exhuast)
hold on
plot(theta, intake)
title('Smoothed Lift Curves')
legend('Exhuast', 'Intake')
ylabel('Lift (mm)')
xlabel('Crank Angle')
