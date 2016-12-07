function [shape] = shape_functions1D(type,gauss)

if strcmp(type,'linear')
    for i=1:gauss.nq
        shape.sf(1,i) = 0.5*(1.0-gauss.pt(i));
        shape.sf(2,i) = 0.5*(1.0+gauss.pt(i));
        shape.dsf(1,i) = -0.5;
        shape.dsf(2,i) = 0.5;
    end
elseif strcmp(type,'quadratic')
    for i=1:gauss.nq
        shape.sf(1,i) = 0.5*(gauss.pt(i)-1.0)*gauss.pt(i);
        shape.sf(2,i) = 1.0-gauss.pt(i)^2;
        shape.sf(3,i) = 0.5*(gauss.pt(i)+1.0)*gauss.pt(i);
        shape.dsf(1,i) = gauss.pt(i)-0.5;
        shape.dsf(2,i) = -2.0*gauss.pt(i);
        shape.dsf(3,i) = gauss.pt(i)+0.5;
    end
else
    error('Mesh type not known')
end
