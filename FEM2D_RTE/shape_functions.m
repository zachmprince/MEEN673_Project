function [shape] = shape_functions(type,gauss)

xnode = [-1  1 1 -1  0 1 0 -1 0; 
         -1 -1 1  1 -1 0 1  0 0];
if strcmp(type,'linear')
    for n=1:4
        k=0;
        for i=1:gauss.nq
            for j=1:gauss.nq
                k=k+1;
                shape.sf(n,k) = 0.25*(1+xnode(1,n)*gauss.pt(i))*(1+xnode(2,n)*gauss.pt(j));
                shape.dsfx(n,k) = 0.25*xnode(1,n)*(1+xnode(2,n)*gauss.pt(j));
                shape.dsfy(n,k) = 0.25*xnode(2,n)*(1+xnode(1,n)*gauss.pt(i));
            end
        end
    end
elseif strcmp(type,'quadratic')
    for n=1:9
        k=0;
        for i=1:gauss.nq
            for j=1:gauss.nq
                k=k+1;
                xi = gauss.pt(i);
                eta = gauss.pt(j);
                xp = xnode(1,n);
                yp = xnode(2,n);
                xi0 = 1.0+xi*xp;
                eta0 = 1.0+eta*yp;
                xi1 = 1.0-xi^2;
                eta1 = 1.0-eta^2;
                xi2 = xp*xi;
                eta2 = yp*eta;
                if n<=4
                    shape.sf(n,k) = 0.25*xi0*eta0*xi2*eta2;
                    shape.dsfx(n,k) = 0.25*xp*eta2*eta0*(1.0+2.0*xi2);
                    shape.dsfy(n,k) = 0.25*yp*xi2*xi0*(1.0+2.0*eta2);
                elseif n==5 || n==7
                    shape.sf(n,k) = 0.5*xi1*eta0*eta2;
                    shape.dsfx(n,k) = -xi*eta2*eta0;
                    shape.dsfy(n,k) = 0.5*xi1*yp*(1.0+2.0*eta2);
                elseif n==6 || n==8
                    shape.sf(n,k) = 0.5*eta1*xi0*xi2;
                    shape.dsfx(n,k) = 0.5*eta1*xp*(1.0+2.0*xi2);
                    shape.dsfy(n,k) = -eta*xi2*xi0;
                else
                    shape.sf(n,k) = xi1*eta1;
                    shape.dsfx(n,k) = -2.0*xi*eta1;
                    shape.dsfy(n,k) = -2.0*eta*xi1;
                end
            end
        end
    end        
else
    error('Mesh type not known')
end

