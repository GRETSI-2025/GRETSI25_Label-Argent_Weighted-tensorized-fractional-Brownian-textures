function s = skewness(x,dim)

if nargin < 2
    dim = find(size(x) ~= 1, 1);
end

x0 = x - mean(x,dim);
std = sqrt(mean(x0.^2,dim));
m3 = mean(x0.^3,dim);
s = m3 ./ std.^3;
% correct the bias ?
%n = size(x,dim);
%s = s .* sqrt((n-1)./n) .* n./(n-2);
end