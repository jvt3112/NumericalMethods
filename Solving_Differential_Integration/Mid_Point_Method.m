clear all;
clc;
L = 30;
m = 68.1;
cd = 0.25;
k = 40;
gamma = 8;
g = 9.81;
% x(0)=0 and v(0)=0
n = 510;
v = zeros(n,1);
vhalf = zeros(n,1);
x = zeros(n,1);
time = zeros(n,1);
time_step = 0.1; % 0.1 seconds 
time(1) = 0;
v(1) = 0;
vhalf(1)=0;
x(1) = 0; %dx/dt = v ; x = v*t
for i=2:n
    time(i) = time(i-1)+time_step;
    cord = 0;
    if(x(i-1)>L)
        cord = (k/m)*(x(i-1)-L) + (gamma/m)*(v(i-1));
    end
    v(i) = v(i-1) + ((g-(sign(v(i-1))*cd*(v(i-1)^2))/m - cord)*time_step)/2;
    x(i) = x(i-1) + v(i-1)*time_step;
    %predictor step
    cord1 = 0;
    if(x(i)>L)
        cord1 = (k/m)*(x(i)-L) + (gamma/m)*(v(i));
    end
    vhalf(i)= (g-(sign(v(i))*cd*(v(i)^2))/m - cord1);
    v(i) = v(i-1) + vhalf(i)*time_step;
    x(i) = x(i-1) + v(i-1)*time_step;
    
end

plot(time,-x);
hold on;
plot(time,v);
legend('x', 'v')
