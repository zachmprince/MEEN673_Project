function [gauss] = quadrature(type)

% Determine appropriate number of points for given element type
if strcmp(type,'linear')
    nq = 2;
elseif strcmp(type,'quadratic')
    nq = 3;
else
    error('Mesh type not known')
end

% Data for gauss points
gauspt=[0.0, 0.0, 0.0, 0.0, 0.0;
%     -0.57735027, 0.57735027, 0.0, 0.0, 0.0;
    -1.0, 1.0, 0.0, 0.0, 0.0;
    -0.77459667, 0.0, 0.77459667, 0.0, 0.0;
    -0.86113631,-0.33998104, 0.33998104, 0.86113631, 0.0;
    -0.90617984,-0.53846931,0.0,0.53846931,0.90617984];

% Data for gauss weights
gauswt=[2.0 0.0 0.0 0.0 0.0 ;
    1.0, 1.0, 0.0, 0.0, 0.0;
    0.55555555, 0.88888888, 0.55555555, 0.0, 0.0;
    0.34785485, 0.65214515, 0.65214515, 0.34785485, 0.0d0;
    0.23692688, 0.47862867, 0.56888888, 0.47862867, 0.23692688];

% Output appropriate parameters
gauss.nq = nq;
gauss.pt = gauspt(nq,1:nq);
gauss.wt = gauswt(nq,1:nq);


