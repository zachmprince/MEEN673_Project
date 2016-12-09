clear;clc;close all

% Mesh paramters
mesh2d.type = 'linear'; % Mesh type
mesh2d.ny = 1; % Number of elements in the y direction
mesh2d.dimx = 1; % Domain width
mesh2d.dimy = 1; % Domain height
time.ntimes = 30; % Number of time steps
time.dt = 0.05;
time.theta = 1.0;

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
param.z = @(x,y) 1.0;
param.k = 0.1;

% Nonlinear Solver Parameters
nonln.type = 'Direct';
nonln.nls = 1;
nonln.sls = 0.0;
nonln.itmax = 100;
nonln.eps = 1e-4;
nonln.beta = 0.5; 
    
% Solve (32X32)
mesh2d.nx = 32; % Number of elements in the x direction
[sol,~] = fem_RTE2D(mesh2d,time,bc,sol,param,nonln);

figure(1)
hold on
plot(sol.x(1,:),sol.E(1,:,end).^0.25,'-',sol.x(1,:),sol.T(1,:,end),'--')
legend('Radiation Temperature','Material Temperature','Location','Best')
xlabel('Distance')
ylabel('Temperature')
axis([0 1 0 1.5])

figure(2)
hold on
plot(sol.x(1,:),sol.E(1,:,end).^0.25,'-')

% Solve (64X64)
mesh2d.nx = 64; % Number of elements in the x direction
[sol,~] = fem_RTE2D(mesh2d,time,bc,sol,param,nonln);

figure(2)
hold on
plot(sol.x(1,:),sol.E(1,:,end).^0.25,'-')

% Solve (128X128)
mesh2d.nx = 128; % Number of elements in the x direction
[sol,~] = fem_RTE2D(mesh2d,time,bc,sol,param,nonln);

figure(2)
hold on
plot(sol.x(1,:),sol.E(1,:,end).^0.25,'-')

% Solve (256X256)
mesh2d.nx = 256; % Number of elements in the x direction
[sol,~] = fem_RTE2D(mesh2d,time,bc,sol,param,nonln);

figure(2)
hold on
plot(sol.x(1,:),sol.E(1,:,end).^0.25,'-')



legend('32X32 Mesh','64X64 Mesh','128X128 Mesh','256X256 Mesh','Location','Best')
xlabel('Distance')
ylabel('Radiation Temperature (E)^1^/^4')
axis([0 1 0 1.5])


