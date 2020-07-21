clear all;
clc;
alpha=0.835; %thermal diffusivity of alluminium 
Lx=40; % length 
Ly=40; %width
dx=2; % space step
dy=2; % space step
nx=Lx/dx+1; % number of intervals
ny=Ly/dy+1;% number of intervals
dt=1; %time_step
t_f=10; %total_time
%boundry conditions
T1=75;  %left wall
T2=50;  %right wall 
T3=0;   %bottom wall
T4=100; %top wall
T=zeros(nx,ny,1); 
tic %starts the time
% Boundary conditions
for i=1:nx
    T(i,1,1)=T1;
    T(i,ny,1)=T2;
end
for j=1:ny
    T(1,j,1)=T3;
    T(nx,j,1)=T4;
end
for t=1:((t_f/dt))
    
    for i=1:nx
        T(i,1,t)=T1;
        T(i,ny,t)=T2;
    end
    for j=1:ny
        T(1,j,t)=T3;
        T(nx,j,t)=T4;
    end
    
    for i=2:(nx-1)
        for j=2:(ny-1)
            T(i,j,t+1)=alpha*dt*((T(i+1,j,t)-2*T(i,j,t)+T(i-1,j,t))/(dx^2)+(T(i,j+1,t)-2*T(i,j,t)+T(i,j-1,t))/(dy^2))+T(i,j,t);
        end
    end
    x=linspace(0,Lx,nx);
    y=linspace(0,Ly,ny);
    surf(y,x,T(:,:,t));   
    axis([0 Ly 0 Lx 0 100]);
    eval(['print -djpeg heat2d_' num2str(t) '.jpeg']);
end
elapsedTime=toc;
fprintf('Took %0.4f seconds \n',elapsedTime);