clear all ; %clears everything 
clc ; %clears command window

n = 11; %number of iterations (0 to 10 min)
Ta = 21; %degree celcius
k = 0.017; % per minute
To = 68; %degree celcius

time = zeros(n,1);
temp = zeros(n,1); %temperature based on euler method
tempexact = zeros(n,1); %temperature based on analytical method

time(1)=0;
temp(1)=To;
tempexact(1)=To;

time_step = 1; % 1 minute 

for i=2:n
    time(i)=time(i-1)+time_step;
    tempexact(i)= Ta + (To-Ta)*exp((-k)*time(i));
    temp(i)=temp(i-1)-(k*(temp(i-1)-Ta)*time_step);
    
end

plot(time,tempexact);
hold on;
plot(time,temp);

legend('Analytical Solution', 'Eulers Method')