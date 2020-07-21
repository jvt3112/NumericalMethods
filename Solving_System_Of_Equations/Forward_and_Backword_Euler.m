clear all;
clc;

dt = 0.1;
n = 500;

% initializinf of arrays
T = zeros(1, n); % time
C1 = zeros(1, n); % array of C1
C2 = zeros(1, n); % array of C2
C3 = zeros(1, n); % array of C3

% initial given values
C1(1) = 1;
C2(1) = 1;
C3(1) = 0;

% Euler's implicit method for solving Stiff ODEs
for i=1:n-1
    T(i+1) = T(i) + dt; 
    %euler's forward method
    C1(i+1) = C1(i)/(1 + (0.013)*dt + 1000*(C3(i))*dt); 
    C2(i+1) = C2(i)/(1 + 2500*(C3(i))*dt); 
    C3(i+1) = (C3(i) - (0.013)*C1(i)*dt)/(1 + 1000*(C1(i))*dt + 2500*(C2(i))*dt);
    
    %euler's backward method
    C1(i+1) = C1(i)/(1 + (0.013)*dt + 1000*(C3(i+1))*dt); 
    C2(i+1) = C2(i)/(1 + 2500*(C3(i+1))*dt); 
    C3(i+1) = (C3(i) - (0.013)*C1(i+1)*dt)/(1 + 1000*(C1(i+1))*dt + 2500*(C2(i+1))*dt);
    
end

C1(n)
C2(n)
C3(n)

% Graph plots
plot(T, C1);
hold on
plot(T, C2);
hold on
plot(T, C3);
grid on
box on
legend('C1', 'C2', 'C3');