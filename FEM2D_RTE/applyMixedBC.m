function [Ke,Fe] = applyMixedBC(Ke,Fe,gauss,shape,elem,bc,ndf)

[npedge,num_edges] = size(elem.bdy);
index=[];
for i=1:num_edges
    if sum(elem.bdy(:,i))>eps
        index = [index i];
    end
end

for ne=1:length(index)
    edge = index(ne);
    for v=1:ndf
        if abs(bc.val1(v,edge))<eps, continue; end
        for n=1:npedge
            if edge==1 || edge==2
                x(n) = elem.y(elem.bdy(n,edge));
            else
                x(n) = elem.x(elem.bdy(n,edge));
            end
        end
        for ni=1:gauss.nq
            gsf = shape.sf(:,ni);
            const = dot(shape.dsf(:,ni),x)*gauss.wt(ni);
            for i=1:npedge
                i0=ndf*(elem.bdy(i,edge)-1)+v;
                Fe(i0) = Fe(i0) + bc.val2(v,edge)/bc.val1(v,edge)*gsf(i)*const;
                for j=1:npedge
                    j0=ndf*(elem.bdy(j,edge)-1)+v;
                    Ke(i0,j0) = Ke(i0,j0) + bc.val0(v,edge)/bc.val1(v,edge)*gsf(i)*gsf(j)*const;
                end
            end
        end
    end
end


