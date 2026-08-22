function [C_eff, P, GH_matrix] = membraneStiffness(ystar, Pj, L, dim, group)
% 
% Based on the paper: "Algorithmic design and effective membrane stiffness of origami and kirigami tessellations and tubes" 
%
% By: Andrew Song, Antoine Moats, Yingchao Peng
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Updated Date: 07/28/26.
%
% Inputs:
% ystar: minimized y coordinate array (3*n by 1 where n is the number of indices)
% Pj: a column cell array of the set of all nodes within each panel 
% L: matrix of the rigidity constraints (L; chi_1) see Eq. (38)
% dim: 2 for planar metamaterials; 3 if otherwise
% group: "translation" or "helical" for group class
%
% Outputs:
% C_eff: effective membrane stiffness tensor in Voigt notation (3x3)
% P: projector matrix used for computing displacement field (nxn)
% GH_matrix: matrix used to determine Guest-Hutchinson modes (NxN) where N = 2I-2(|B1|+|B2|+1)

textOutput = 0;
tol = 10^(-10); % for rounding numerics

% initialize vectors
e1 = [1; 0; 0];
e2 = [0; 1; 0];
e3 = [0; 0; 1]; % out of plane vector
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

% compute mu_ij 2D or 3D
if dim == 2

    for j = 1:lenJ
        % compute dummy sums for omega_jstar
        num = zeros(n,1); % numerator for omega_jstar; see Eq. (49, 51)
        den = 0; % denominator for omega_jstar; Eq. (49, 51)
        for i = 1:length(Pj{j}) 
            k = Pj{j}(i); % vertex k
            num = num + chi_ij{k,j}'*cross(e3, y_ijstar{k,j}); 
            den = den + norm(cross(e3, y_ijstar{k,j}))^2;
        end
        omega_jstar{j} = num/den;

        % compute mu_ij
        for i = 1:length(Pj{j})
            k = Pj{j}(i); % vertex k
            mu_ij{k,j} = chi_ij{k,j} - cross(e3, y_ijstar{k,j})*omega_jstar{j}';
        end

    end
else % 3D case
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
end

for j = 1:lenJ
    for i = 1:length(Pj{j})
        k = Pj{j}(i); % vertex k
        Astar = Astar + mu_ij{k,j}'*mu_ij{k,j}; % in this case, Astar is mu'*mu
    end
end

% clean Astar
Astar = (Astar+Astar'); % symmetricize and normalize for m = 1/2 u dot Astar u (correction for Astar = 2mu'*mu)
Astar(abs(Astar)<tol) = 0; % round off numerics

N_G = null(L); % null space of (L; chi_1) sized 3I by NG

% clean N_G'*Astar*N_G
GH_matrix = N_G'*Astar*N_G;
GH_matrix = (GH_matrix+GH_matrix')/2; % symmetricize
GH_matrix(abs(GH_matrix)<tol) = 0; % round off numerics

% check for G-H mode special case
GH_det = det(GH_matrix);
if GH_det < tol % N_G'*Astar*N_G is singular
    %disp("Guest-Hutchinson mode found.")
    %fprintf("%d Guest-Hutchinson modes found.\n", size(null(GH_matrix),2))
end

% compute projector and stiffness matrix using pseudoinverse
P = eye(n) - N_G*pinv(GH_matrix)*N_G'*Astar; 
Kstar = Astar*P;
% clean Kstar
Kstar = (Kstar+Kstar')/2; % symmetricize
Kstar(abs(Kstar)<tol) = 0; % round off numerics

Ystar = constructY(ystar, group);

C_eff = Ystar'*Kstar*Ystar; % Voigt notation representation of the effective membrane stiffness tensor

% clean C_eff
C_eff = (C_eff+C_eff')/2; % symmetricize
C_eff(abs(C_eff)<tol) = 0; % round off numerics

% Mechanism check
[V,D] = eig(C_eff);
zero_idx = abs(diag(D)) < tol;
num_mech = sum(zero_idx);
if textOutput
    if num_mech > 0
        fprintf("%d DoF mechanism \n", num_mech)
        disp("Strain space modes: ")
        disp(V(:, zero_idx))
    else, disp("Bistability demonstrated.")
    end
end

end