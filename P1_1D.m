% Clear existing constants, close all figures, clear console
clear; 
close all; 
clf;

% Parameters
D  = 0.5;         % diffusion coefficient
T = 1;            % simulation time duration
L = 5;            % length of 1D domain

% Spatio-temporal domain discretization
dx = 0.01;
x = -L:dx:L;
Nx = length(x);

dt = 0.00009; 
Nt = T/dt;
t = linspace(0,T,Nt);

% Initial condition
u       = zeros(length(x),length(t));
u(:,1)  = 25*exp(-10*x.^2);   % initial temperature distribution in degrees Celsius

% Add a buried landmine
landmine_position = 0;        % position of the landmine
landmine_radius = 0.2;        % radius of the landmine
landmine_depth = 0.5;         % depth of the landmine
u(x>=landmine_position-landmine_radius & x<=landmine_position+landmine_radius,1) = 75;  % set temperature of landmine to 75 degrees Celsius

% Compute the 1D discrete Laplacian for interior points
A = (D*dt/dx^2)*lap1d(Nx);
    
% Modify the 1D discrete Laplacian to account for boundary points
% Homogeneous Neumann BC (du/dx = 0 at x = -L and x = L)
A(1,1) = -1;
A(Nx,Nx) = -1;

outevery = 10;
for i = 1:Nt
    if (mod(i,outevery)==0)
        plot(x,u(:,i),'-m','linewidth',1);
        ylim([0 100]);       % set y-axis limits to 0 and 100 degrees Celsius
        title(sprintf('Temperature Distribution over Time (t = %8.6f seconds)',(i-1)*dt),'fontweight','normal');
        xlabel('Position');
        ylabel('Temperature (Celsius)');
        set(gca,'fontname','CMU Serif'); box on; grid on;
        set(gca,'FontSize',20);
        set(gcf,'color','w');
        hold off;
        currFrame = getframe(gcf);
    end
    u(:,i+1) = u(:,i) + A*u(:,i);
end

% lap1d
%
% Form the (scaled) 1D Laplacian for Dirichlet boundary conditions
% on a node-centered grid
%
% Input: N -- number of grid points (no bdy pts)
%
% Output: L -- N x N sparse matrix for discrete Laplacian

function L = lap1d(N)
    e=ones(N,1);
    L = spdiags([e -2*e e],-1:1,N,N);
end