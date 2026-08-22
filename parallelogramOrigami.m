function [x, Pj, B, B_k] = parallelogramOrigami(n, vec)
% parallelogramOrigami: a function to construct parallelogram panel origami
% nodal positions and panel labeling with n by n unit cells given 2*n edge vectors
%
% By: Andrew Song
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Updated Date: 08/17/26.
%
% Inputs:
% n: (even) number of unit cells along each axis
% vec: 3 by n^2 sized array of vector columns parameterizing edge lengths
% ex) n = 2:
% v_1 + ... + v_n = l1R
% v_n+1 + ... + v_2n = l2R
%
% Outputs:
% x: a column vector of nodal positions, sized 3*(n+1)^2
% Pj: a column cell array of the set of all nodal indices within each panel (the jth panel
% corresponds to the jth cell)
% B: a matrix of pairs of nodal indices describing edge connectivity, sized (2n(n+1), 2)
% B_k: a matrix of pairs of nodal indices related by periodicity, sized (|B1|+|B2|, 2)

N = n + 1;

%% Construct nodal positions x
x = zeros(3*N^2,1);
for r = 0:n
    % coordinate of leftmost node per row
    % 1st node
    if r == 0, xLeft = [0;0;0];
    else
        % use leftmost node of previous row
        prevLeft = nodeLabel(n,r-1,0);
        xLeft = x(3*prevLeft-2 : 3*prevLeft);
        for k = 0:n-1
            if mod(r,n) == k
                if k == 0, xLeft = xLeft + vec(:,2*n);
                else, xLeft = xLeft + vec(:,n+k);
                end
            end
        end
    end

    % construct row: left -> right
    xRow = zeros(3,N);
    xRow(:,1) = xLeft;
    for c = 1:n
        for k = 0:n-1
            if mod(c,n) == k
                if k == 0, xRow(:,c+1) = xRow(:,c) + vec(:,n);
                else, xRow(:,c+1) = xRow(:,c) + vec(:,k); end
            end
        end
    end

    for c = 0:n
        node = nodeLabel(n,r,c);
        x(3*node-2 : 3*node) = xRow(:,c+1);

    end
end

%% Define panel labeling Pj using a snake pattern
J = n^2;
Pj = cell(J,1);

for r = 0:n-1
    for c = 0:n-1
        
        n1 = nodeLabel(n,r,c);
        n2 = nodeLabel(n,r,c+1);
        n3 = nodeLabel(n,r+1,c+1);
        n4 = nodeLabel(n,r+1,c);
        
        % determine panel number using snake pattern        
        if mod(r,2) == 0, j = r*n + c + 1; % left -> right
        else, j = (r+1)*n - c; % left <- right
        end
        
        Pj{j} = sort([n1 n2 n3 n4]);
    end
end

%% Construct edge connectivity matrix B
B = zeros(2*n*N,2);

for i = 1:N^2-1 % for snake pattern
    B(i,:) = [i i+1];
end

k = N^2 - 1;
for r = 0:n-1
    if mod(r,2) == 0 % even rows
        for c = 0:n-1
            node1 = nodeLabel(n,r,c);
            node2 = nodeLabel(n,r+1,c);

            k = k + 1;
            B(k,:) = [node1 node2];
        end
    else % odd rows
        for c = 1:n
            node1 = nodeLabel(n,r,c);
            node2 = nodeLabel(n,r+1,c);

            k = k + 1;
            B(k,:) = [node1 node2];
        end
    end
end

%% Construct matrix of periodicity-related index pairs B_k
B_k = zeros(2*n + 1, 2); % # of periodicity pairs: (n+1) left/right + n bottom/top
k = 1;

% left-right boundary pairs
for r = 0:n
    leftNode  = nodeLabel(n,r,0);
    rightNode = nodeLabel(n,r,n);
    B_k(k,:) = [leftNode rightNode];
    k = k + 1;
end

% bottom-top boundary pairs
for c = 0:n-1 % exclude rightmost column to avoid redundant corner constraint
    bottomNode = nodeLabel(n,0,c);
    topNode    = nodeLabel(n,n,c);
    B_k(k,:) = [bottomNode topNode];
    k = k + 1;
 end

%% Helper function to access node index based on row/column #
function node = nodeLabel(n,r,c)

    N = n + 1;
    if mod(r,2) == 0
        % even row: left -> right
        node = r*N + c + 1;
    else
        % odd row: left <- right
        node = (r+1)*N - c;
    end

end

end