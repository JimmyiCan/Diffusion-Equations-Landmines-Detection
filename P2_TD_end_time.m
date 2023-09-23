%2D Heat Equation.
clear; close all; clc
n = 10;  % grid has 10 points per dimension (overlapping)
x = linspace(-n/2,n/2,n); % x axis
dx = x(2)-x(1); % delta x
y = linspace(-n/2,0,n); % y axix
dy = dx; % delta y
TOL = 1e-6; % tollerence   
T = zeros(n); % Define the temperature matrix
T(1,1:n) = exp(-0.1*x.^2); %Bottom
T(n,1:n) = 0;  %Top
% T(1:n,1) = 0.001*exp(-0.1*(y+n/2).^2);  %LEFT
% T(1:n,n) = 0.001*exp(-0.1*(y+n/2).^2);  %RIGHT
dt = dx^2/4; % delta t
error = 1; 
k = 0;
    while error > TOL 
        k = k+1;
          Told = T;
          for i = 2:n-1
              for j = 2:n-1 % 2d heat equation
                      T(i,j) = dt*((Told(i+1,j)-2*Told(i,j)+Told(i-1,j))/dx^2 ... 
                      + (Told(i,j+1)-2*Told(i,j)+Told(i,j-1))/dy^2) ...
                      + Told(i,j);
              end
          end
          error = max(max(abs(Told-T)));
    end
    
T = 100*T;

subplot(2,1,1)
contour(x,y,T)
title('Temperature (Steady State)','FontSize',18)
xlabel('x','FontSize',14)
ylabel('y','FontSize',14)
colorbar
subplot(2,1,2)
pcolor(x,y,T)
shading interp
title('Temperature (Steady State)','FontSize',18)
xlabel('x','FontSize',14)
ylabel('y','FontSize',14),colorbar

