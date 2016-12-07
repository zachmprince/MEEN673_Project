function [sol,mesh,time,gauss,shape] = generate_mesh_RTE(mesh,time,sol)

% Determine number of nodes and nodes per element
if strcmp(mesh.type,'linear')
    mesh.Nx = mesh.nx + 1;
    mesh.Ny = mesh.ny + 1;
    mesh.npe = 4;
    mesh.npedge = 2;
elseif strcmp(mesh.type,'quadratic')
    mesh.Nx = mesh.nx*2 + 1;
    mesh.Ny = mesh.ny*2 + 1;
    mesh.npe = 9;
    mesh.npedge = 3;
else
    error('Mesh type not known')
end

% Compute gauss quadrature parameters
[gauss] = quadrature(mesh.type);

% Compute shape function values at quadrature points
[shape] = shape_functions(mesh.type,gauss);
[shape.oneD] = shape_functions1D(mesh.type,gauss);


% Allocate arrays
num_elem = mesh.nx*mesh.ny;
sol.x = zeros(mesh.Ny,mesh.Nx);
sol.y = zeros(mesh.Ny,mesh.Nx);
sol.t = time.dt*(0:time.ntimes);
sol.E = zeros(mesh.Ny,mesh.Nx,time.ntimes+1);
sol.T = zeros(mesh.Ny,mesh.Nx,time.ntimes+1);

% Generate arrays for x, y, and initial u accounting for Dirichlet BCs
dx = mesh.dimx/(mesh.Nx-1);
dy = mesh.dimy/(mesh.Ny-1);
for i=1:mesh.Ny
    for j=1:mesh.Nx
        sol.x(i,j) = (j-1)*dx;
        sol.y(i,j) = (i-1)*dy;
        sol.E(i,j,1) = sol.E0;
        sol.T(i,j,1) = sol.T0;
    end
end

% Generate connectivity array:
% Row # is the element #
% Column # is the local node #
% Entries are the corresponding global node #s
mesh.elem = zeros(num_elem,mesh.npe);
n=0;
for i=1:mesh.nx
    for j=1:mesh.ny
        n=n+1;
        if n==1
            mesh.elem(n,1) = 1;
        elseif j>1
            mesh.elem(n,1) = mesh.elem(n-1,4);
        else
            mesh.elem(n,1) = mesh.elem(n-mesh.ny,2);
        end
           
        if strcmp(mesh.type,'linear')
            mesh.elem(n,2) = mesh.elem(n,1) + mesh.Ny;
            mesh.elem(n,3) = mesh.elem(n,2) + 1;
            mesh.elem(n,4) = mesh.elem(n,1) + 1;
        elseif strcmp(mesh.type,'quadratic')
            mesh.elem(n,2) = mesh.elem(n,1) + mesh.Ny*2;
            mesh.elem(n,3) = mesh.elem(n,2) + 2;
            mesh.elem(n,4) = mesh.elem(n,1) + 2;
            mesh.elem(n,5) = mesh.elem(n,1) + mesh.Ny;
            mesh.elem(n,6) = mesh.elem(n,2) + 1;
            mesh.elem(n,7) = mesh.elem(n,5) + 2;
            mesh.elem(n,8) = mesh.elem(n,1) + 1;
            mesh.elem(n,9) = mesh.elem(n,5) + 1;
        else
            error('Mesh type not known')
        end
    end
end

% Generate boundary connectivity array:
% Row # is the element #
% Column # is the local edge #
% Depth # is the linear node #
% Value is local node #
if strcmp(mesh.type,'linear')
    mesh.bdy = zeros(num_elem,4,2);
    n=0;
    for i=1:mesh.nx
    for j=1:mesh.ny
        n=n+1;
        if j==1
            mesh.bdy(n,3,1) = 1;
            mesh.bdy(n,3,2) = 2;
        end
        if j==mesh.ny 
            mesh.bdy(n,4,1) = 4;
            mesh.bdy(n,4,2) = 3;
        end
        if i==1
            mesh.bdy(n,1,1) = 1;
            mesh.bdy(n,1,2) = 4;
        end
        if i==mesh.nx 
            mesh.bdy(n,2,1) = 2;
            mesh.bdy(n,2,2) = 3;
        end    
    end
    end
elseif strcmp(mesh.type,'quadratic')
    mesh.bdy = zeros(num_elem,4,3);
    n=0;
    for i=1:mesh.nx
    for j=1:mesh.ny
        n=n+1;
        if j==1
            mesh.bdy(n,3,1) = 1;
            mesh.bdy(n,3,2) = 5;
            mesh.bdy(n,3,3) = 2;
        end
        if j==mesh.ny 
            mesh.bdy(n,4,1) = 4;
            mesh.bdy(n,4,2) = 7;
            mesh.bdy(n,4,3) = 3;
        end
        if i==1
            mesh.bdy(n,1,1) = 1;
            mesh.bdy(n,1,2) = 8;
            mesh.bdy(n,1,3) = 4;
        end
        if i==mesh.nx 
            mesh.bdy(n,2,1) = 2;
            mesh.bdy(n,2,2) = 6;
            mesh.bdy(n,2,3) = 3;
        end    
    end
    end
end


% Time constants
time.a1 = time.theta*time.dt;
time.a2 = (1.0-time.theta)*time.dt;

    
    