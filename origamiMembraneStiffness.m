function [C_eff, P, GH_matrix] = origamiMembraneStiffness(ystar, Pj, L, X)
% 
% Based on the paper: "Algorithmic design and effective membrane stiffness of origami and kirigami tessellations and tubes" 
% By: Andrew Song, Antoine Moats, Yingchao Peng
% for use in implementing the dynamical system in the 3D origami case discussed in the paper:
% "Dynamical system for the evolution of a mechanism motion in origami and kirigami metamaterials"
% By: Andrew Song, Yingchao Peng, Ian Tobasco, and Paul Plucinsky
%
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Updated Date: 09/25/26.
%
% Inputs:
% ystar: minimized y coordinate array (3*n by 1 where n is the number of indices)
% Pj: a column cell array of the set of all nodes within each panel 
% L: matrix of the rigidity constraints (L; chi_1) see Eq. (38)
% X: a matrix such that Ex = Xe where E is the strain matrix [e1 e3/2; e3/2 d1] and e is the strain vector in Voigt notation
%
% Outputs:
% C_eff: effective membrane stiffness tensor in Voigt notation (3x3)
% P: projector matrix used for computing displacement field (nxn)
% GH_matrix: matrix used to determine Guest-Hutchinson modes (NxN) where N = 2I-2(|B1|+|B2|+1)

% initialize vectors
n = length(ystar); % number of vertices * 3
lenJ = length(Pj); % number of panels
y_ijstar = cell(n/3, lenJ); % stores information of R_jstar*x_ijstar or approximately y_istar - avg(ystar)_j at an energy minimum
[chi, ~] = construct_chi(ystar); % chi{i}*y = y_i
chi_ij = cell(n/3, lenJ); % chi_ij{i,j} stores information of chi_i - avg(chi)_j
omega_jstar = cell(lenJ, 1); % intermediate tensor used for mu computing
mu_ij = cell(3, n); % intermediate tensor used for Astar computing
Astar = zeros(n, n); % intermediate tensor used for Kstar computing

% construct y_ijstar and chi_ij
for j = 1:lenJ

    num_Pj = length(Pj{j}); % number of vertices in the j-th panel
    chi_avg = (1/num_Pj)*calcMatrixSum_y(chi, Pj{j});

    for i = 1:length(Pj{j})
        k = Pj{j}(i); % vertex k
        chi_ij{k,j} = chi{k} - chi_avg; % j dependency encoded in chi_avg
        y_ijstar{k,j} = chi_ij{k,j}*ystar;
    end

end

% compute mu_ij^(3D)
% compute dummy sums for omega_jstar
for j = 1:lenJ

    left = zeros(3); % tensor meant to be inverted for omega_jstar; see Eq. (52, 54)
    right = zeros(3, n); % right tensor in omega_jstar; see Eq. (52, 54)
    for i = 1:length(Pj{j})
        k = Pj{j}(i); % vertex k
        y_ijstarcross = [0                   -y_ijstar{k,j}(3)    y_ijstar{k,j}(2);
            y_ijstar{k,j}(3)     0                   -y_ijstar{k,j}(1);
            -y_ijstar{k,j}(2)    y_ijstar{k,j}(1)     0];
        left = left + y_ijstarcross'*y_ijstarcross;
        right = right + y_ijstarcross'*chi_ij{k,j};
    end
    omega_jstar{j} = left\right;

    % compute mu_ij
    for i = 1:length(Pj{j})
        k = Pj{j}(i); % vertex k
        y_ijstarcross = [0                   -y_ijstar{k,j}(3)    y_ijstar{k,j}(2);
            y_ijstar{k,j}(3)     0                   -y_ijstar{k,j}(1);
            -y_ijstar{k,j}(2)    y_ijstar{k,j}(1)     0];
        mu_ij{k,j} = chi_ij{k,j} - y_ijstarcross*omega_jstar{j};
    end
end

for j = 1:lenJ
    for i = 1:length(Pj{j})
        k = Pj{j}(i); % vertex k
        Astar = Astar + mu_ij{k,j}'*mu_ij{k,j}; % in this case, Astar is mu'*mu
    end
end

% clean Astar
Astar = (Astar+Astar'); % symmetricize and normalize for m = 1/2 u dot Astar u (correction for Astar = 2mu'*mu)

N_G = null(L); % null space of (L; chi_1) sized 3I by NG

% clean N_G'*Astar*N_G
GH_matrix = N_G'*Astar*N_G;
GH_matrix = (GH_matrix+GH_matrix')/2; % symmetricize

% compute projector and stiffness matrix using pseudoinverse
P = eye(n) - N_G*(GH_matrix\(N_G'*Astar)); 
Kstar = Astar*P;
% clean Kstar
Kstar = (Kstar+Kstar')/2; % symmetricize

C_eff = X'*Kstar*X;

% clean C_eff
C_eff = (C_eff+C_eff')/2; % symmetricize

end