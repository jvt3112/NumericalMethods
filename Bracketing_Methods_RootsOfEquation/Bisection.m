clear all;
clc;

k1 = 50000; % g/s^2
g = 9.81; % m/s^2
m = 90; % g
k2 = 40; % g/s^2*m^1/2
h = 0.45;% m 

% decide the bracket 

x_left = 0;
x_right = 5;

mid = (x_left + x_right)/2;

fxmid = ((2*k2*(mid^(5/2)))/5) + ((1/2)*k1*(mid^2)) - m*g*mid - m*g*h;
while(abs(fxmid)>1e-03)
    fxl = ((2*k2*(x_left^(5/2)))/5) + ((1/2)*k1*(x_left^2)) - m*g*x_left - m*g*h;
    fxmid = ((2*k2*(mid^(5/2)))/5) + ((1/2)*k1*(mid^2)) - m*g*mid - m*g*h;
    if((fxl*fxmid)==0)
        break;
    end
    if((fxl*fxmid)<0) 
        x_right = mid;
        mid = (x_left + x_right)/2;
    else
        x_left = mid;
        mid = (x_left + x_right)/2;
    end
end





