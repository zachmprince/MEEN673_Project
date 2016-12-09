function [Ke,Fe] = element_mat2D_RTE(gauss,shape,elem,param,newton,time,bc)

npe = length(elem.x);

Ke = zeros(2*npe,2*npe); 
% Ke0 = zeros(2*npe,2*npe);
Ce = zeros(2*npe,2*npe); 
Fe = zeros(2*npe,1);

if newton
    Te = zeros(2*npe,2*npe); 
end

k=0;
for ni=1:gauss.nq
for nj=1:gauss.nq
    k=k+1;
    % Compute Jacobian
    Jac = [shape.dsfx(:,k)'; shape.dsfy(:,k)']*[elem.x elem.y];
    Jd = det(Jac);
    
    % Transform shape functions
    gsf = shape.sf(:,k)';
    gdsf = Jac\[shape.dsfx(:,k)'; shape.dsfy(:,k)'];
    const = Jd*gauss.wt(ni)*gauss.wt(nj);

    x = dot(elem.x',gsf);
    y = dot(elem.y',gsf);
    z = param.z(x,y);
    T = dot(elem.T',gsf);
    E = dot(elem.E',gsf);
    Ex = dot(elem.E',gdsf(1,:));
%     Tp = dot(elem.Tp',gsf);

    if newton
        E = dot(elem.E',gsf);
        Ex = dot(elem.E',gdsf(1,:));
        Ey = dot(elem.E',gdsf(2,:));
        Tx = dot(elem.T',gdsf(1,:));
        Ty = dot(elem.T',gdsf(2,:));
    end

    % Assemble coefficient arrays
    ii=1;
    for i=1:npe
        jj=1;
        for j=1:npe
            scc = gsf(i)*gsf(j)*const;
            sxx = gdsf(1,i)*gdsf(1,j)*const;
            syy = gdsf(2,i)*gdsf(2,j)*const;
            
            Ke(ii,jj) = Ke(ii,jj) + T^3/(3*z^3)*(sxx+syy) + z^3/T^3*scc;
            Ke(ii,jj+1) = Ke(ii,jj+1) - z^3*param.a*scc;
            Ke(ii+1,jj) = Ke(ii+1,jj) - param.c*z^3/T^3*scc;
            Ke(ii+1,jj+1) = Ke(ii+1,jj+1) + param.k*T^(5/2)*(sxx+syy) + param.c*z^3*param.a*scc;
            
%             Ke011(i,j) = Ke011(i,j) + K11_func(Tp,scc,sxx,syy,z);
%             Ke012(i,j) = Ke012(i,j) + K12_func(scc,z);
%             Ke021(i,j) = Ke021(i,j) + K21_func(Tp,scc,z);
%             Ke022(i,j) = Ke022(i,j) + K22_func(Tp,scc,sxx,syy,z);
            
            Ce(ii,jj) = Ce(ii,jj) + 1/param.c*scc;
            Ce(ii+1,jj+1) = Ce(ii+1,jj+1) + param.Cv*scc;
            
            if newton
                sxc = gdsf(1,i)*gsf(j)*const;
                syc = gdsf(2,i)*gsf(j)*const;
                Te(ii,jj+1) = Te(ii,jj+1) + T^2/z^3*(sxc*Ex+syc*Ey) - 3*z^3/T^4*E*scc;
                Te(ii+1,jj+1) = Te(ii+1,jj+1) + 5/2*T^(3/2)*(sxc*Tx+syc*Ty) + 3*param.c*z^3/T^4*E*scc;
            end
            jj=2*j+1;
        end
        ii=2*i+1;
    end
end
end

if sum(sum(elem.bdy)) > eps
    [Ke,Fe] = applyMixedBC(Ke,Fe,gauss,shape.oneD,elem,bc,2);
%     [Ke0,~] = applyMixedBC(Ke,Fe,gauss,shape.oneD,elem,bc,2);
end

sol = zeros(2*npe,1); solp = sol;
ii=1;
for i=1:npe
    sol(ii)   = elem.E(i);
    sol(ii+1) = elem.T(i);
    solp(ii)   = elem.Ep(i);
    solp(ii+1) = elem.Tp(i);
    ii=2*i+1;
end

% Fe = (time.a2*Ke0+Ce)*solp + (time.a1+time.a2)*Fe;
Fe = (Ce)*solp + (time.a1+time.a2)*Fe;
Ke = time.a1*Ke + Ce;

if newton
    Fe = Fe - Ke*sol;
    Ke = Ke + time.a1*Te;
end


