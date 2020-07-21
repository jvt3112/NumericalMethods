y = [0 1.8 2 4 4 6 4 3.6 3.4 2.8 0];
n = 10; %number of sub intervals
a = 0; %lower limit value f(a) = f(0) = y(1) = 0
lower_value = 0;
b = 20; % higher limit value f(b) = f(20) = y(11) = 0
higher_value = 0;
h = (b-a)/n; 
so=0;se=0;
for k=2:1:n
    if rem(k,2)==1
       so=so+y(k);%sum of odd terms
     else
       se=se+y(k); %sum of even terms
    end
end
so
se
% Formula:  (h/3)*[(y1+yn)+2*(y3+y5+..odd term)+4*(y2+y4+y6+...even terms)]
answer=h/3*(lower_value + higher_value + 4*se+2*so);
fprintf('\n The value of integration is %f',answer);