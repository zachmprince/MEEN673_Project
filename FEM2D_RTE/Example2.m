clear;clc;close all

% Mesh paramters
mesh2d.type = 'linear'; % Mesh type
mesh2d.nx = 16; % Number of elements in the x direction
mesh2d.ny = 16; % Number of elements in the y direction
mesh2d.dimx = 1; % Domain width
mesh2d.dimy = 1; % Domain height

% Time paramters
time.type = 'Parabolic'; % None, Parabolic, or Hyperbolic
time.gamma = 1.0; % Number of elements in the y direction

% Boundary and Initial Conditions
%  bc.val0*val + bc.val1*D*grad(val) = bc.val2;
bc.val0 = [1/4 1/4 0 0; ...% [left right bottom top]
           0 0 0 0];   % [left right bottom top]
bc.val1 = [-1/2 1/2 1 1; ...% [left right bottom top]
           1 1 1 1];   % [left right bottom top]
bc.val2 = [1 0 0 0; ...% [left right bottom top]
           0 0 0 0];   % [left right bottom top]
sol.E0 = 1e-5;
sol.T0 = sol.E0^0.25;

% Equation Parameters
param.Cv = 1.0;
param.c = 1.0;
param.a = 1.0;
param.z = @(x,y) 1.0 + 9.0*(x>=1/3)*(x<=2/3)*(y>=1/3)*(y<=2/3);
param.k = 0.1;

% Nonlinear Solver Parameters
nonln.type = 'Direct';
nonln.nls = 1;
nonln.sls = 0.0;
nonln.itmax = 100;
nonln.eps = 1e-4;
nonln.beta = 0.5; 
    
% Solve (BE)
time.ntimes = 60; % Number of time steps
time.dt = 0.05; % Number of elements in the x direction
time.theta = 1.0; % Number of elements in the x direction
[sol,~] = fem_RTE2D(mesh2d,time,bc,sol,param,nonln);

for i=1:12
    figure(i)
    hold on
    % surf(sol.x,sol.y,sol.E(:,:,end).^0.25)
    contour(sol.x,sol.y,sol.T(:,:,i*5),'ShowText','on')
end

