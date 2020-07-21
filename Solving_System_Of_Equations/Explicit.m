clear all;
clc;

Cp = 0.2174;
ro = 2.7;
k = 0.49;
delta_x = 1; % total length 10 cm so division 0 2 4 6 8 10
delta_t = 0.5; % total time 9 sec so division 1 2 3 4 5 6 7 8 9 
total_length=10;
total_time = 10; 
xdiv = total_length/delta_x + 1;
tdiv = total_time/delta_t;
% gamma = alpha*(delta_t)/(delta_x)^2
alpha = k/(Cp*ro) ;
gamma = alpha*(delta_t)/(delta_x)^2 ;

T = zeros(tdiv,xdiv);
T(:,1) = 100; 
T(:,xdiv) = 50;

for i=2:tdiv-1
    for j=2:xdiv-1
        T(i,j) = T(i-1,j) + gamma*(T(i-1,j+1)- 2*(T(i-1,j)) + T(i-1,j-1));
    end
end
T
plot(T(2:xdiv-1,:)')



