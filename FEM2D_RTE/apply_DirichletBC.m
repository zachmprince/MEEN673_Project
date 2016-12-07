function [K, F] = apply_DirichletBC(K,F,u,mesh,bc,newton)

n = 0;
for j=1:mesh.Nx
for i=1:mesh.Ny
for k=1:2
    n=n+1;
    if j==1 % Left boundary
        if abs(bc.val1(k,1))<eps 
            K(n,:) = zeros(size(K(n,:)));
            K(n,n) = 1;
            if newton
                F(n) = bc.val2(k,1)/bc.val0(k,1) - u(n);
            else
                F(n) = bc.val2(k,1)/bc.val0(k,1);
            end
        end
    end    
    if j==mesh.Nx % Right boundary
        if abs(bc.val1(k,2))<eps 
            K(n,:) = zeros(size(K(n,:)));
            K(n,n) = 1;
            if newton
                F(n) = bc.val2(k,2)/bc.val0(k,2) - u(n);
            else
                F(n) = bc.val2(k,2)/bc.val0(k,2);
            end
        end
    end    
    if i==1 % Bottom boundary
        if abs(bc.val1(k,3))<eps 
            K(n,:) = zeros(size(K(n,:)));
            K(n,n) = 1;
            if newton
                F(n) = bc.val2(k,3)/bc.val0(k,3) - u(n);
            else
                F(n) = bc.val2(k,3)/bc.val0(k,3);
            end
        end
    end    
    if i==mesh.Ny % Top boundary
        if abs(bc.val1(k,4))<eps 
            K(n,:) = zeros(size(K(n,:)));
            K(n,n) = 1;
            if newton
                F(n) = bc.val2(k,4)/bc.val0(k,4) - u(n);
            else
                F(n) = bc.val2(k,4)/bc.val0(k,4);
            end
        end
    end
end
end
end
