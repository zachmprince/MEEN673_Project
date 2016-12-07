function [K,F] = assemble2D_RTE(mesh,time,sol,param,gauss,shape,nt,newton,bc)

E_curr = sol.E(:,:,nt+1);
E_prev = sol.E(:,:,nt);
T_curr = sol.T(:,:,nt+1);
T_prev = sol.T(:,:,nt);

% Allocate global arrays
num_nodes = mesh.Nx*mesh.Ny;
num_elem = mesh.nx*mesh.ny;
K = sparse([],[],[],2*num_nodes,2*num_nodes,2*mesh.npe*num_elem);
F = zeros(2*num_nodes,1);

% Assemble arrays
elem.x  = zeros(mesh.npe,1);
elem.y  = zeros(mesh.npe,1);
elem.E  = zeros(mesh.npe,1);
elem.Ep = zeros(mesh.npe,1);
elem.T  = zeros(mesh.npe,1);
elem.Tp = zeros(mesh.npe,1);
elem.bdy = zeros(mesh.npedge,4);
for n=1:num_elem
    % Determine elemental values of x, y, and u
    for i=1:mesh.npe
        elem.x(i)  = sol.x(mesh.elem(n,i));
        elem.y(i)  = sol.y(mesh.elem(n,i));
        elem.E(i)  = E_curr(mesh.elem(n,i));
        elem.Ep(i) = E_prev(mesh.elem(n,i));
        elem.T(i)  = T_curr(mesh.elem(n,i));
        elem.Tp(i) = T_prev(mesh.elem(n,i));
    end
    for i=1:mesh.npedge
        for j=1:4
            elem.bdy(i,j) = mesh.bdy(n,j,i);
        end
    end
    
    % Compute element arrays
    [Ke,Fe] = element_mat2D_RTE(gauss,shape,elem,param,newton,time,bc);
    
    % Assemble element arrays to global arrays
    for i=1:mesh.npe
        for ni=1:2
            F((mesh.elem(n,i)-1)*2+ni) = F((mesh.elem(n,i)-1)*2+ni) + Fe((i-1)*2+ni);
            for j=1:mesh.npe
                for nj=1:2
                    K((mesh.elem(n,i)-1)*2+ni,(mesh.elem(n,j)-1)*2+nj) = K((mesh.elem(n,i)-1)*2+ni,(mesh.elem(n,j)-1)*2+nj) + Ke((i-1)*2+ni,(j-1)*2+nj);
                end
            end
        end
    end
end
        









