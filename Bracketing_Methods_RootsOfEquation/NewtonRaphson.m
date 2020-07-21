clear all;
clc;

k1 = 50000; % g/s^2
g = 9.81; % m/s^2
m = 90; % g
k2 = 40; % g/s^2*m^1/2
h = 0.45;% m 

xg = 0.1;
fxg = ((2*k2*(xg^(5/2)))/5) + ((1/2)*k1*(xg^2)) - m*g*xg - m*g*h;

while(abs(fxg)>0.001)
    
    fxg = ((2*k2*(xg^(5/2)))/5) + ((1/2)*k1*(xg^2)) - m*g*xg - m*g*h;
    fdashxg = k2*(xg^(3/2)) + k1*xg -m*g;
    xg = xg - (fxg/fdashxg);
end
