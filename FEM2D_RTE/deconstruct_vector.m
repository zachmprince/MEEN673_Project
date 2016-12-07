function [varargout] = deconstruct_vector(u,varargin)

ndf = numel(varargin);
num_nodes = numel(varargin{1});
varargout = cell(size(varargin));
for i=1:ndf
    varargout{i} = zeros(size(varargin{i}));
end

for n=1:num_nodes
    for v=1:ndf
        nn = ndf*(n-1)+v;
        varargout{v}(n) = u(nn);
    end
end
