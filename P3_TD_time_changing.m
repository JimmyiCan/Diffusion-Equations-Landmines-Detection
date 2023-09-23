%2D Heat Equation.
clear; close all; clc
n = 10;  % grid has 10 points per dimension (overlapping)
x = linspace(-n/2,n/2,n); % x axis
dx = x(2)-x(1); % delta x
y = linspace(-n/2,0,n); % y axix
dy = dx; % delta y
TOL = 1e-6; % tollerence   
T = zeros(n); % Define the temperature matrix


T(1,1:n) = exp(-0.1*1*x.^2); %Bottom
T(n,1:n) = 0;  %Top
% T(1:n,1) = 0.001*exp(-0.1*(y+n/2).^2);  %LEFT
% T(1:n,n) = 0.001*exp(-0.1*(y+n/2).^2);  %RIGHT
D=1;
dt = dx^2/8*D; % delta t
error = 1; 
k = 0;


figure

% Define the resolution of the video
set(gcf,'Position',[50 50 7680 4320]) % 8K
% set(gcf,'Position',[50 50 3840 3160]) % 4K
% set(gcf,'Position',[50 50 2560 1440]) % 2K
% set(gcf,'Position',[50 50 1920 1080]) % 1080p, full HD, most common
% set(gcf,'Position',[50 50 1280 720]) % 720p, HD
% set(gcf,'Position',[50 50 854 480]) % 480p


v = VideoWriter('Single Landmine.mp4','MPEG-4');
v.Quality = 100; % range = [0,100]
open(v);

 
grid on ; 


while error > TOL 
      k = k+1;
      Told = T;
      for i = 2:n-1
          for j = 2:n-1 % 2d heat equation
              T(i,j) = D*dt*((Told(i+1,j)-2*Told(i,j)+Told(i-1,j))/dx^2 ... 
                      + (Told(i,j+1)-2*Told(i,j)+Told(i,j-1))/dy^2) ...
                      + Told(i,j);
          end
      end
      error = max(max(abs(Told-T)));
      figure(1)
      pcolor(x,y,T);
      title('Temperature','FontSize',18)
      xlabel('x','FontSize',14)
      ylabel('y','FontSize',14)
      shading interp;
      getframe(gcf);


      frame = getframe(gcf);
      writeVideo(v,frame);

end

close(v)

T = 100*T;
subplot(2,1,1)
contour(x,y,T)
title('Temperature (Steady State)','FontSize',18)
xlabel('x','FontSize',14)
ylabel('y','FontSize',14)
colorbar

subplot(2,1,2)
pcolor(x,y,T);
shading interp
title('Temperature (Steady State)','FontSize',18)
xlabel('x','FontSize',14)
ylabel('y','FontSize',14)
colorbar



