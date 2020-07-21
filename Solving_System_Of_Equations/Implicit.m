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
solution_matrix = zeros(xdiv-2,1);
T(2,2) = gamma*100;
T(2,xdiv-1) = gamma*50;
solution_matrix(1,1) = gamma*100;
solution_matrix(xdiv-2,1) = gamma*50;
solution_matrix
for i=2:tdiv-1
    A = zeros(xdiv-2,xdiv-2);
    for j=2:xdiv-3
        A(j,j-1) = -gamma;
        A(j,j) = 1+2*gamma;
        A(j,j+1) = -gamma;
    end
    A(1,1) = 1+2*gamma;
    A(1,2) = -gamma;
    A(xdiv-2,xdiv-3) = -gamma;
    A(xdiv-2,xdiv-2) = 1+2*gamma;
    A
    ans = TDMA(A,solution_matrix,gamma);
    for k=2:xdiv-1
        T(i+1,k) = ans(k-1,1);
        solution_matrix(k-1,1) = T(i+1,k);
    end
    solution_matrix;
end
T
plot(T(2:xdiv-1,:)')

function y = TDMA(A,B,gamma)
B(1,1) = B(1,1)+gamma*100;
B(end,1)=B(end,1) + gamma*50;
f = B;
n = length(f);
v = zeros(n,1);   
y = v;
a = diag(A);
b = zeros(n,1);
b(2:end) = diag(A,-1);
c = zeros(n,1);
c(1:end-1) = diag(A,+1);
w = a(1);
y(1) = f(1)/w;
for i=2:n
    v(i-1) = c(i-1)/w;
    w = a(i) - b(i)*v(i-1);
    y(i) = ( f(i) - b(i)*y(i-1) )/w;
end
for j=n-1:-1:1
   y(j) = y(j) - v(j)*y(j+1);
end
end

