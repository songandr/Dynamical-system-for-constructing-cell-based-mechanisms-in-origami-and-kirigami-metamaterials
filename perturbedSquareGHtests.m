% This script is designed to run test cases for the perturbed 2x2 square lattice
% structure discussed in the paper:
%
% "Dynamical system for constructing cell-based mechanisms in periodic metamaterials"
% by Andrew Song, Yingchao Peng, Ian Tobasco, and Paul Plucinsky
%
% By: Andrew Song
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California
%
% Updated Date: 09/23/26.
%
% This test is interested in the Guest-Hutchinson modes of
% a planar 2x2 symmetric square lattice structure with a perturbed center.

close all
clear
clc
cla

% Initialize all reference configurations

% symmetric state
x_sym = [0 0;
        1 0;
        2 0;
        2 1;
        1 1;
        0 1;
        0 2;
        1 2;
        2 2];

% generate perturbations along each axis
n = 401;
deltax_center_1 = linspace(-1,1,n);
deltax_center_2 = linspace(-1,1,n);

x_perturbed = cell(n);
for i = 1:n
    for j = 1:n
        perturbation = zeros(9,2);
        deltax_center = [deltax_center_1(i), deltax_center_2(j)];
        perturbation(5,:) = deltax_center;
        x_perturbed{i,j} = x_sym + perturbation;
        x_perturbed{i,j} = reshape(x_perturbed{i,j}',[],1);
    end
end

B = [1 2;
    2 3;
    3 4;
    4 5;
    5 6;
    6 7;
    7 8;
    8 9;
    2 5;
    1 6;
    5 8;
    4 9];

B_k = [1 3;
    6 4;
    7 9;
    1 7;
    2 8];

% Construct null space matrix N
L_0 = periodicityMatrix(B_k, 9, 2);
N = null(L_0); % null space matrix, size 2I by N where N = 2I-2(|B1|+|B2|+1)

% Construct GH matrices N^TKN and find their eigenvalues
lambda_2 = zeros(n); % 2nd smallest eigenvalues
lambda_m = zeros(n); % minimum eigenvalues
for i = 1:n
    for j = 1:n
        K_0 = barStiffness(x_perturbed{i,j},B);
        GH_0 = N'*K_0*N;
        GH_sym = (GH_0+GH_0')/2; % symmetricize
        GH_eigs = eig(GH_sym);
        lambda_m(i,j) = GH_eigs(1);
        lambda_2(i,j) = GH_eigs(2);
    end
end

%% Plotting

% Initialize
tol = 1e-14;
[X,Y] = meshgrid(deltax_center_1, deltax_center_2);

% Construct colormap
blue  = [0.00 0.25 0.80];
white = [0.98 0.98 0.98];
red   = [0.80 0.00 0.00];
nblue = 128;
nred = 127;
bluemap = [linspace(blue(1), white(1), nblue)', ...
           linspace(blue(2), white(2), nblue)', ...
           linspace(blue(3), white(3), nblue)'];
redmap = [linspace(white(1), red(1), nred)', ...
          linspace(white(2), red(2), nred)', ...
          linspace(white(3), red(3), nred)'];
cmap = [bluemap; redmap(2:end,:)]; % blue to white to red, no white duplicate
ncolors = size(cmap,1);

% Plot lambda_min
lambda_m = lambda_m';
lambda_m_max = max(lambda_m(:));  % maximum value
nonzero = lambda_m > tol;
nonzero_lambda_m = lambda_m(nonzero); % nonzero values
nonzero_lambda_m_min = min(nonzero_lambda_m); % smallest nonzero value
log_min = log10(nonzero_lambda_m_min);
log_max = log10(lambda_m_max);
log_white = (log_min + log_max)/2; % white is the mid point

color_index = ones(size(lambda_m));
log_lambda_m = log10(lambda_m(nonzero));
color_index(nonzero) = ...
    1 + (ncolors-1) * ...
    (log_lambda_m - log_min) / (log_max-log_min);
color_index(nonzero) = max(1, min(ncolors,round(color_index(nonzero)))); % white is the midpoint
color_index(~nonzero) = 1; % values < tol are blue

figure
surf(X, Y, zeros(size(X)), color_index, ...
    'EdgeColor', 'none');
view(2)
axis equal tight
set(gca, 'XTickLabel', {}, 'YTickLabel', {}, 'ZTickLabel', {});
colormap(cmap)
clim([1 ncolors])
%xlabel('\delta_1')
%ylabel('\delta_2')
%title('Minimum Eigenvalues of Perturbed 2x2 Square Lattices')

cb = colorbar;
%cb.Label.String = '\lambda_m(N^T K(x)N)';

powers = -5:-1;
tick_values = 10.^powers;
tick_positions = zeros(1,length(powers));
tick_positions(1) = 1;
tick_positions(2:end) = 1 + (ncolors-1) * (log10(tick_values(2:5))-log_min) / (log_max-log_min);
%cb.Ticks = tick_positions;
%cb.TickLabels = arrayfun(@(p) sprintf('10^{%d}', p), powers, 'UniformOutput', false);
cb.TickLabels = [];

% Plot lambda_2
lambda_2 = lambda_2';
lambda_2_max = max(lambda_2(:));  % maximum value
nonzero2 = lambda_2 > tol;
nonzero_lambda_2 = lambda_2(nonzero); % nonzero values
nonzero_lambda_2_min = min(nonzero_lambda_2); % smallest nonzero value
log_min2 = log10(nonzero_lambda_2_min);
log_max2 = log10(lambda_2_max);
log_white2 = (log_min2 + log_max2)/2; % white is the mid point

color_index2 = ones(size(lambda_2));
log_lambda_2 = log10(lambda_2(nonzero2));
color_index2(nonzero2) = ...
    1 + (ncolors-1) * ...
    (log_lambda_2 - log_min2) / (log_max2-log_min2);
color_index2(nonzero2) = max(1, min(ncolors,round(color_index2(nonzero2)))); % white is the midpoint
color_index2(~nonzero2) = 1; % values < tol are blue

figure
surf(X, Y, zeros(size(X)), color_index2, ...
    'EdgeColor', 'none');
view(2)
axis equal tight
set(gca, 'XTickLabel', {}, 'YTickLabel', {}, 'ZTickLabel', {});
colormap(cmap)
clim([1 ncolors])
%xlabel('\delta_1')
%ylabel('\delta_2')
%title('Second Eigenvalues of Perturbed 2x2 Square Lattices')

cb2 = colorbar;
%cb2.Label.String = '\lambda_2(N^T K(x)N)';

tick_positions2 = zeros(1,length(powers));
tick_positions2(1) = 1;
tick_positions2(2:end) = 1 + (ncolors-1) * (log10(tick_values(2:5))-log_min2) / (log_max2-log_min2);
%cb2.Ticks = tick_positions2;
%cb2.TickLabels = arrayfun(@(p) sprintf('10^{%d}', p), powers, 'UniformOutput', false);
cb2.TickLabels = [];

%% Functions
function K = barStiffness(y, B)

% construct chi_ij
chi = construct_chi2D(y); % chi{i}*y = y_i

% compute stiffness matrix K
K = zeros(length(y)); % 2I by 2I matrix
for b = 1:size(B,1) % goes down the list of index pairs
    i = B(b,1); % first index
    j = B(b,2); % second index
    chi_ij = chi{i} - chi{j};
    y_ij = chi_ij*y;
    n = y_ij/norm(y_ij);
    K = K + chi_ij'*(n*n')*chi_ij;
end

% clean K
K = (K+K')/2; % symmetricize
end