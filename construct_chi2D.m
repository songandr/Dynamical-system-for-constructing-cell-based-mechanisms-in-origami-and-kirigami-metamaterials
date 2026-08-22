function chi = construct_chi2D(x)
% 
% constructs the map chi s.t. chi_i*x = x_i
%
% Inputs
% x: x coordinate 2-D column array (2*n by 1 where n is the number of vertices)
% 
% Outputs
% chi: cell array where i-th entry is the 2 by 2*n mapping matrix
% mapping x in R^(2n) to x_i in R^2

% Initialize relevant parameters
n = length(x); % number of vertices * 2
chi = cell(n/2,1);

for k = 1:(n/2)
    chi_k = zeros(2, n);
    chi_k(1:2, 2*k-1:2*k) = eye(2);
    chi{k} = chi_k;
end