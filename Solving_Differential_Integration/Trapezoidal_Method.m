clc;
clear all;
x = [0 0.02 0.04 0.05 0.06 0.07 0.1]; % in meters
r = [0.002 0.00135 0.00134 0.0016 0.00158 0.00142 0.002]; %in meters
n = 6; %number of intervals
f=@(x)1/(x^4);

mew = 0.005; % N-s/m^2
Q = 10*(10^(-6)); % m^3/s

%using trapezoidal rule
ans = 0;
for k=2:7
   ans = ans + ((x(k)-x(k-1))*(f(r(k))+f(r(k-1)))/2);
end
(ans*(Q)*(-8)*(mew))/pi
