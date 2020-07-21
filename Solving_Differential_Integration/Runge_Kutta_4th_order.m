clear all;
clc;
% x(0)=0 and v(0)=0
n = 51;
v = zeros(n,1);
x = zeros(n,1);
time = zeros(n,1);
h = 1; % 0.1 seconds 
time(1) = 0;
v(1) = 0;
x(1) = 0; %dx/dt = v ; x = v*t

% 4th order runge-kutta methos 
% y(i+1) = y(i) + (k1 + 2*k2 + 2*k3 + k4)/6
% k1 = h*f(x(i),y(i))
% k2 = h*f(x(i)+h/2,y(i)+k1/2)
% k3 = h*f(x(i)+h/2,y(i)+k2/2)
% k4 = h*f(x(i)+h ,y(i) + k3)

for i=2:n
    time(i) = time(i-1)+h;
    k_1 = runge(x(i-1),v(i-1));
    k_2 = runge(x(i-1)+0.5*h,v(i-1)+0.5*k_1*h);
    k_3 = runge((x(i-1)+0.5*h),(v(i-1)+0.5*k_2*h));
    k_4 = runge((x(i-1)+h),(v(i-1)+k_3*h));
    v(i) = v(i-1) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    x(i) = x(i-1) + v(i-1)*h;
    
end
plot(time,-x);
hold on;
plot(time,v);
legend('x', 'v')

function [output]=runge(pos,vel)
    L = 30;
    m = 68.1;
    cd = 0.25;
    k = 40;
    gamma = 8;
    g = 9.81;
    cord = 0;
    if(pos>L)
        cord = (k/m)*(pos-L) + (gamma/m)*(vel);
    end
    output = g-sign(vel)*cd*((vel)^2)/m - cord;
end