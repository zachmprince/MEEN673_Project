function u = construct_vector(varargin)

ndf = numel(varargin);
num_nodes = numel(varargin{1});
u = zeros(ndf*num_nodes,1);

for n=1:num_nodes
    for v=1:ndf
        nn = ndf*(n-1)+v;
        u(nn) = varargin{v}(n);
    end
end

