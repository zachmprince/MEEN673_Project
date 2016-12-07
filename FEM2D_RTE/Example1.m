clear;clc;close all

% Mesh paramters
mesh2d.type = 'linear'; % Mesh type
mesh2d.nx = 8; % Number of elements in the x direction
mesh2d.ny = 3; % Number of elements in the y direction
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
bc.val2 = [-1 0 0 0; ...% [left right bottom top]
           0 0 0 0];   % [left right bottom top]
sol.E0 = 1e-5;
sol.T0 = sol.E0^0.25;

% Equation Parameters
param.Cv = 1.0;
param.c = 1.0;
param.a = 1.0;
param.z = 1.0;
param.k = 1.0;

% Nonlinear Solver Parameters
nonln.type = 'Newton';
nonln.nls = 1;
nonln.sls = 0.0;
nonln.itmax = 10;
nonln.eps = 1e-3;
nonln.beta = 0; 
    
% Solve (BE)
time.ntimes = 100; % Number of time steps
time.dt = 0.001; % Number of elements in the x direction
time.theta = 1.0; % Number of elements in the x direction
[sol,~] = fem_RTE2D(mesh2d,time,bc,sol,param,nonln);

figure(1)
surf(sol.x,sol.y,sol.E(:,:,end))

figure(2)
surf(sol.x,sol.y,sol.T(:,:,end))
