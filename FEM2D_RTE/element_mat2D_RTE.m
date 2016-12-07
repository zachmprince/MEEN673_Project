function [Ke,Fe] = element_mat2D_RTE(gauss,shape,elem,param,newton,time,bc)

npe = length(elem.x);

Ke = zeros(2*npe,2*npe); 
Ke0 = zeros(2*npe,2*npe);
Te = zeros(2*npe,2*npe); 
Ce = zeros(2*npe,2*npe); 
Fe = zeros(2*npe,1);

Ke11 = zeros(npe,npe); Ke12 = Ke11; Ke21 = Ke11; Ke22 = Ke11; 
Ke011 = Ke11; Ke012 = Ke11; Ke021 = Ke11; Ke022 = Ke11; 
Te12 = Ke11; Te22 = Ke11;
Ce11 = Ke11; Ce22 = Ke11;
Fe1 = zeros(npe,1); Fe2 = Fe1;

K11_func = @(T,scc,sxx,syy) T^3/(3*param.z^3)*(sxx+syy) + param.z^3/T^3*scc;
K12_func = @(scc) -param.z^3*param.a*scc;
K21_func = @(T,scc) -param.c*param.z^3/T^3*scc;
K22_func = @(T,scc,sxx,syy) param.k*abs(T)^(5/2)*(sxx+syy) + param.c*param.z^3*param.a*scc;

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

    T = dot(elem.T',gsf);
    Tp = dot(elem.Tp',gsf);

    if newton
        E = dot(elem.E',gsf);
        Ex = dot(elem.E',gdsf(1,:));
        Ey = dot(elem.E',gdsf(2,:));
        Tx = dot(elem.T',gdsf(1,:));
        Ty = dot(elem.T',gdsf(2,:));
    end

    % Assemble coefficient arrays
    for i=1:npe
        for j=1:npe
            scc = gsf(i)*gsf(j)*const;
            sxx = gdsf(1,i)*gdsf(1,j)*const;
            syy = gdsf(2,i)*gdsf(2,j)*const;
            
            Ke11(i,j) = Ke11(i,j) + K11_func(T,scc,sxx,syy);
            Ke12(i,j) = Ke12(i,j) + K12_func(scc);
            Ke21(i,j) = Ke21(i,j) + K21_func(T,scc);
            Ke22(i,j) = Ke22(i,j) + K22_func(T,scc,sxx,syy);
            
            Ke011(i,j) = Ke011(i,j) + K11_func(Tp,scc,sxx,syy);
            Ke012(i,j) = Ke012(i,j) + K12_func(scc);
            Ke021(i,j) = Ke021(i,j) + K21_func(Tp,scc);
            Ke022(i,j) = Ke022(i,j) + K22_func(Tp,scc,sxx,syy);
            
            Ce11(i,j) = Ce11(i,j) + 1/param.c*scc;
            Ce22(i,j) = Ce22(i,j) + param.Cv*scc;
            
            if newton
                sxc = gdsf(1,i)*gsf(j)*const;
                syc = gdsf(2,i)*gsf(j)*const;
                Te12(i,j) = Te12(i,j) + T^2/param.z^3*(sxc*Ex+syc*Ey) - 3*param.z^3/T^4*E*scc;
                Te22(i,j) = Te22(i,j) + 5/2*T^(3/2)*(sxc*Tx+syc*Ty) + 3*param.c*param.z^3/T^3*scc;
            end
        end
    end
end
end

% Rearrangment
ii=1;
for i=1:npe
    Fe(ii) = Fe1(i);
    Fe(ii+1) = Fe2(i);
    jj=1;
    for j=1:npe
        Ke(ii,jj)     = Ke11(i,j);
        Ke(ii,jj+1)   = Ke12(i,j);
        Ke(ii+1,jj)   = Ke21(i,j);
        Ke(ii+1,jj+1) = Ke22(i,j);
        Ke0(ii,jj)     = Ke011(i,j);
        Ke0(ii,jj+1)   = Ke012(i,j);
        Ke0(ii+1,jj)   = Ke021(i,j);
        Ke0(ii+1,jj+1) = Ke022(i,j);
        Ce(ii,jj)     = Ce11(i,j);
        Ce(ii+1,jj+1) = Ce22(i,j);
        Te(ii,jj+1)   = Te12(i,j);
        Te(ii+1,jj+1) = Te22(i,j);
        jj=2*j+1;
    end
    ii=2*i+1;
end

if sum(sum(elem.bdy)) > eps
    [Ke,Fe] = applyMixedBC(Ke,Fe,gauss,shape.oneD,elem,bc,2);
    [Ke0,~] = applyMixedBC(Ke,Fe,gauss,shape.oneD,elem,bc,2);
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

Fe = (time.a2*Ke0+Ce)*solp + (time.a1+time.a2)*Fe;
Ke = time.a1*Ke + Ce;

if newton
    Fe = Fe - Ke*sol;
    Ke = Ke + time.a1*Te;
end


