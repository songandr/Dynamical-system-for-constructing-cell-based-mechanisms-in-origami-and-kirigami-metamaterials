function L = periodicityMatrix(B,num_nodes,dim)
    
% B: ordered node index pairs, size m x 2
% num_nodes: total number of nodes
% dim: spatial dimension, 2 or 3
%
% Each row B(k,:) = [i i'] imposes x_i' - x_i = 0.
% The final dim rows anchor node 1:
%   [I_dim  0  ...  0]

m = size(B,1);

% Initialize constraint matrix
L = zeros(dim*(m+1), dim*num_nodes);

% Periodicity constraints
for k = 1:m

    i = B(k,1);
    iprime = B(k,2);

    % Block corresponding to node i
    L(dim*(k-1)+1 : dim*k, dim*(i-1)+1 : dim*i) = -eye(dim);

    % Block corresponding to node i'
    L(dim*(k-1)+1 : dim*k, dim*(iprime-1)+1 : dim*iprime) = eye(dim);
end

% Anchor node 1
L(dim*m+1 : dim*(m+1), 1:dim) = eye(dim);

end