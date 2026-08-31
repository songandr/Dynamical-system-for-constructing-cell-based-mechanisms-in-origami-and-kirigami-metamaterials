% This script is designed as a test case for the MATLAB function 
% mechanismODE and its integration which is based on the paper:
% 
% "Dynamical system for the evolution of a mechanism motion in origami
% and kirigami metamaterials" by Paul Plucinsky and Ian Tobasco.
%
% By: Andrew Song
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Updated Date: 08/31/26.
%
% These tests are interested in the mechanism motion of:
% 1. a 4 bar quad unit cell (1 DoF mechanism)
% 2. a 2D 2x2 quad lattice (1 DoF mechanism..?)
% 3. a rigid rotating squares pattern (1 DoF mechanism)
% 4. a 2x2 Kagome lattice (twisted: 1 DoF, otherwise: GH)
% 5. a 2x2 quad Miura-ori unit cell (1 DoF mechanism)
% 6. a nxn parallelogram unit cell
% 7. a 2x2 quad eggbox unit cell (1 DoF mechanism)
% 8. a "Morph" origami; see Geometric Mechanics of Origami Patterns Exhibiting Poisson’s Ratio Switch
%                       by Breaking Mountain and Valley Assignment (PP Pratapa, Ke Liu, GH Paulino) 2019

close all
clear
clc
cla

% User-defined inputs
tol = 10^(-15);
plotMechanism = 1; % for mechanism animation
tessellate = 0; % 2x2 unit cell tessellation for plotting
snapshotPlots = 0;
eigPlot = 1;
energyPlot = 1;
strainPlot = 0;
residualCheck = 1;
reverse = 0; % default setting; some tests flip initial eigenvector sign for desired mechanism motion
test = 1; % 1: 1x1 square lattice (simplest mechanism test)
          % 2: 2x2 quad lattice (dynamics can't run unless perturbed center)
          % 3: 2x2 rotating squares (kirigami mechanism test)
          % 4: 1x1 twisted Kagome lattice (GH mode at most open state)
          % 5: 2x2 parallelogram origami (rigid panel energy test)
          % 6: nxn parallelogram origami (general quad origami)
          % 7: eggbox origami
          % 8: "Morph" origami
    
% x is the IC of y, a column vector of length 2I holding nodal positions
if test == 1 % 1x1 square lattice
    reverse = 1;
    %theta = 0.01*(pi/2); % pi/2 corresponds to the square state & 0 the flat state
    %T = 126.35; % PSD threshold crossed at t = 126.354561
    reverse = 0;
    theta = 1*(pi/2);
    T = 1;
    x = [0 0;
         1 0;
         1+cos(theta) sin(theta);
         cos(theta) sin(theta)];
    x = reshape(x',[],1);
    x(abs(x)<tol) = 0;
    
    % B is an array of node indices connected by bars
    B = [1 2;
         3 4;
         1 4;
         2 3]; % for a 4 bar unit cell
    
    B_k = [1 2;
           4 3;
           1 4];

    % L_0 is the corresponding periodic conditions in 2D, stacked with chi_1
    L_0 = periodicityMatrix(B_k, 4, 2);

elseif test == 2 % 2x2 square lattice
    epsilon = zeros(9, 2); % symmetric case
    epsilon(5,:) = [0.25 0.4]; % perturbed case
    %epsilon = [-0.2 0.3];

    % line shifted case 1
    %{
    epsilon(2,:) = [0 0.2]; 
    epsilon(5,:) = [0 0.2];
    epsilon(8,:) = [0 0.2];
    %}
    % line shifted case 2
    %{
    epsilon(4,:) = [0.2 0];
    epsilon(5,:) = [0.2 0];
    epsilon(6,:) = [0.2 0];
    %}
    % generic case
    %{
    epsilon(2,:) = [0.1 0.2];
    epsilon(6,:) = [0.15 0.1];
    epsilon(5,:) = epsilon(2,:) + epsilon(6,:);
    epsilon(3,:) = [0.2 -0.3];
    epsilon(4,:) = epsilon(3,:) + epsilon(5,:) - epsilon(2,:);
    epsilon(7,:) = [-0.25 0.1];
    epsilon(8,:) = epsilon(7,:) + epsilon(5,:) - epsilon(6,:);
    epsilon(9,:) = epsilon(4,:) + epsilon(8,:) - epsilon(5,:);
    %}
    x = [0 0;
         1 0;
         2 0;
         2 1;
         1 1;
         0 1;
         0 2;
         1 2;
         2 2];
    x = x + epsilon;
    x = reshape(x',[],1);
    
    T = 1.3;
    T = 1;
    reverse = 0;

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

    % L_0 is the corresponding periodic conditions in 2D, stacked with chi_1
    L_0 = periodicityMatrix(B_k, 9, 2);

elseif test == 3 % rotating squares
    reverse = 0;
    xi = 0.05*pi/4; % 0: fully closed, 1: fully open
    T = 0.511; % 1 GH mode at t=0.511669
    s = 0.5 / (cos(xi)+sin(xi)); % side length given |l_1| = 1
    x1 = [0; 0];
    x2 = s*[cos(xi); -sin(xi)];
    x3 = s*[sin(xi)+cos(xi); cos(xi)-sin(xi)];
    x4 = s*[2*sin(xi)+cos(xi); -sin(xi)];
    x5 = s*[2*sin(xi)+2*cos(xi); 0];
    x6 = s*[sin(xi)+2*cos(xi); cos(xi)];
    x7 = s*[2*sin(xi)+2*cos(xi); 2*cos(xi)];
    x8 = s*[2*sin(xi)+cos(xi); 2*cos(xi)+sin(xi)];
    x9 = s*[sin(xi)+cos(xi); cos(xi)+sin(xi)];
    x10 = s*[cos(xi); 2*cos(xi)+sin(xi)];
    x11 = s*[0; 2*cos(xi)];
    x12 = s*[sin(xi); cos(xi)];
    x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12];

    B = [1 2;
         2 3;
         3 4;
         4 5;
         5 6;
         6 7;
         7 8;
         8 9;
         9 10;
         10 11;
         11 12;
         1 12;
         3 12;
         3 6;
         6 9;
         9 12;
         1 3;
         3 5;
         7 9;
         9 11]; % for rigid panels

    B_k = [1 5;
           11 7;
           2 10;
           4 8];
    
    L_0 = periodicityMatrix(B_k, 12, 2);

elseif test == 4 % twisted Kagome lattice
    % twisted test 1
    %{
    x1 = [0; 0];
    x2 = [1; 0];
    x3 = [2; 0];
    x4 = [3; 0];
    x5 = [4; 0];
    x6 = [0.5; sqrt(3)/2];
    x7 = [2.5; sqrt(3)/2];
    x8 = [4.5; sqrt(3)/2];
    x9 = [1; sqrt(3)];
    x10 = [2; sqrt(3)];
    x11 = [3; sqrt(3)];
    x12 = [4; sqrt(3)];
    x13 = [5; sqrt(3)];
    x14 = [1.5; 3*sqrt(3)/2];
    x15 = [3.5; 3*sqrt(3)/2];
    x16 = [5.5; 3*sqrt(3)/2];
    x17 = [2; 4*sqrt(3)/2];
    x18 = [3; 4*sqrt(3)/2];
    x19 = [4; 4*sqrt(3)/2];
    x20 = [5; 4*sqrt(3)/2];
    x21 = [6; 4*sqrt(3)/2];

    x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12; x13; x14; x15; x16; x17; x18; x19; x20; x21];
    B = [1 2;
         2 3;
         3 4;
         4 5;
         9 10;
         10 11;
         11 12;
         12 13;
         17 18;
         18 19;
         19 20;
         20 21;
         1 6;
         6 9;
         9 14;
         14 17;
         3 7;
         7 11;
         11 15;
         15 19;
         5 8;
         8 13;
         13 16;
         16 21;
         2 6;
         4 7;
         7 10;
         10 14;
         8 12;
         12 15;
         15 18;
         16 20];
    L_0 = [                -eye(2)   zeros(2, 2*3)   eye(2)  zeros(2, 2*16);
                           -eye(2)   zeros(2, 2*15)  eye(2)  zeros(2, 2*4);
           zeros(2, 2*8)   -eye(2)   zeros(2, 2*3)   eye(2)  zeros(2, 2*8);
           zeros(2, 2*2)   -eye(2)   zeros(2, 2*15)  eye(2)  zeros(2, 2*2);
           zeros(2, 2*4)   -eye(2)   zeros(2, 2*15)  eye(2);
           zeros(2, 2*5)   -eye(2)   zeros(2)        eye(2)  zeros(2, 2*13);
           zeros(2, 2*13)  -eye(2)   zeros(2)        eye(2)  zeros(2, 2*5);
           zeros(2)        -eye(2)   zeros(2, 2*15)  eye(2)  zeros(2, 2*3);
           zeros(2, 2*3)   -eye(2)   zeros(2, 2*15)  eye(2)  zeros(2);
                                                     eye(2)  zeros(2, 2*20)];
    %}
    
    % twisted test 2
    s = 1;
    theta = pi/2+0.99*pi/6;
    T = 1.4; % 1 GH mode at t=1.401481
    x1 = [0;0];
    x2 = x1+[s;0];
    x3 = x2+s*[cos(theta);-sin(theta)];
    x4 = x3+s*[cos(pi/3);sin(pi/3)];
    x5 = x4+s*[-cos(pi/3+theta); sin(pi/3+theta)];
    x6 = x5-[s;0];
    x7 = x5-(x3-x1);
    x8 = x1+s*[cos(pi/3);sin(pi/3)];
    x = [x1; x2; x3; x4; x5; x6; x7; x8];
   
    B = [1 2;
         2 3;
         3 4;
         4 5;
         5 6;
         6 7;
         7 8;
         1 8;
         2 8;
         6 8];
    
    B_k = [1 3;
           8 4;
           7 5;
           1 6;
           2 5];
    
    L_0 = periodicityMatrix(B_k, 8, 2);
    

elseif test == 5 % 2x2 Miura origami
    n = 2;
    shear = 0;
    flattening = 1;
    vec = [[1; 0.2; 0], [1; -0.2; 0], [0; 1; 1], [0; 1; -1]]/2;
    [x,Pj,B,B_k] = parallelogramOrigami(n, vec); % standard Miura-ori axial shape change
    if flattening, reverse = 0; T = 0.31; % T = 0.346 until flat state for flattening
    else, T = 0.4; reverse = 1;
    end

    if shear % Miura-ori shear (gamma = 0.5) shape change
    x = [0;0;0;
        0.636828212556466;0.301327900039711;1.51674513879982e-05;
        0.902368927062183;0;0;
        1.29418604239371;0.223710694618478;0.217593713193951;
        1.02860453621252;0.525027541972545;0.217608511651365;
        0.391817115331531;0.223710694618478;0.217593713193951;
        0.783611624891224;0.447462290995219;0;
        1.42043983744769;0.748790191034930;0;
        1.68598055195341;0.447462290995219;0];
    T = 3.4; % flat state (PSD threshold crossed at t = 0.535073)
    end
    
    L = periodicityMatrix(B_k, 9, 3);

elseif test == 6 % n by n parallelogram origami
    n = 6;
    T = 0.5;
    reverse = 1;
    l1R = [1; 0; 0];
    l2R = [0; 1; 0];
    
    % Randomized ICs
    rng(3, "twister");
    r_theta_p = normrnd(pi/4, pi/12, [n, 1]); % 1st n/2 for l1, 2nd n/2 for l2
    r_theta_m = normrnd(-pi/4, pi/12, [n-2, 1]); % 1st n/2-1 for l1, 2nd n/2-1 for l2
    r_theta_p = r_theta_p + [zeros(n/2,1); pi/2*ones(n/2,1)];
    r_theta_m = r_theta_m + [zeros(n/2-1,1); pi/2*ones(n/2-1,1)];

    vec = zeros(3, 2*n);
    for i=1:n-1
        if mod(i,2) == 1
            theta = r_theta_p((i+1)/2);
        else
            theta = r_theta_m(i/2);
        end
        vec(1:2,i) = [cos(theta); sin(theta)]/n/(sqrt(2)/2); % avg 1/n
    end
    vec(:,n) = l1R - sum(vec(:, 1:n-1),2); % enforce l1R
    
    for i=1:n-1
        if mod(i,2) == 1
            theta = r_theta_p(n/2+(i+1)/2);
        else
            theta = r_theta_m(n/2-1+i/2);
        end
        vec(1:2,n+i) = [cos(theta); sin(theta)]/n/(sqrt(2)/2); % avg 1/n
    end
    vec(:,2*n) = l2R - sum(vec(:, n+1:2*n-1),2); % enforce l2R

    r_z_mountain = normrnd(0.5, 0.2, [1, n/2]);
    r_z_valley = normrnd(-0.5, 0.2, [1, n/2-1]);

    for i=1:n/2
        idx_m = n + 2*i - 1; % Mountain edge
        vec(3,idx_m) = r_z_mountain(i);
    
        idx_v = n + 2*i; % Valley edge
        if i < n/2
            vec(3,idx_v) = r_z_valley(i);
        else
            % Last vector closes the lattice in z.
            vec(3,idx_v) = -sum(vec(3,n+1:idx_v-1));
        end
    end

    [x,Pj,B,B_k] = parallelogramOrigami(n, vec); % standard Miura-ori axial shape change
    
    L = periodicityMatrix(B_k,(n+1)^2,3);

elseif test == 7 % eggbox
    flattening = 1;
    vec = [[1; 0; 0.5], [1; 0; -0.5], [0; 1; 0.5], [0; 1; -0.5]]/2;
    [x,Pj,B,B_k] = parallelogramOrigami(2, vec);
    
    if flattening, T = 0.346; reverse = 1; % until flat state for flattening
    else, T = 3.58; reverse = 0;
    end
    
    L = periodicityMatrix(B_k, 9, 3);

elseif test == 8 % Morph
    T = 1;
    reverse = 1; % 1: transition; 0: maintain M or E mode
    %hybrid_state = "EEEEMMMM";

    % Miura: alpha + beta = pi; eggbox: beta = alpha;
    alpha = pi/3;
    beta = 4*pi/18;
    %beta = pi-alpha;
 
    % Morph configurational space is fully described by:
    % phi(0 <= psi <= psi_max = 2*beta < pi)
    % psi(0 < phi_min = alpha - beta <= phi <= phi_max = alpha + beta < pi)
    phi = 8*pi/18;
    psi = acos(cos(2*alpha) + 2*(cos(beta) - cos(alpha)*cos(phi))^2/sin(phi)^2); 
    disp("psi_max = "+rad2deg(2*beta))
    disp("psi = "+rad2deg(psi))

    a = 1;
    c = 1;
    b = a*abs(cos(alpha)/cos(beta)); % orthorhombic Morph cell condition
    % Reference Fig. 2a of PP Pratapa, K. Liu, GH Paulino (2019)
    % L = sqrt(a^2+b^2-2*a*b*cos(phi));
    % W = 2*c*sin(psi/2);

    v1 = [a*sqrt(cos(psi/2)^2 - cos(alpha)^2)/cos(psi/2); 0; a*cos(alpha)/cos(psi/2)]; % O1O2
    v2 = [b*sqrt(cos(psi/2)^2 - cos(beta)^2)/cos(psi/2); 0; -b*cos(beta)/cos(psi/2)]; % O2O3
    v3 = [0; c*sin(psi/2); c*cos(psi/2)]; % O1O4
    v4 = [0; c*sin(psi/2); -c*cos(psi/2)]; % O4O7

    x1 = [0; 0; 0];
    x2 = x1 + v1;
    x3 = x2 + v2;
    x6 = x1 + v3; % O4
    x5 = x2 + v3; % O5
    x4 = x3 + v3; % O6
    x7 = x6 + v4;
    x8 = x5 + v4;
    x9 = x4 + v4;    
    x = [x1; x2; x3; x4; x5; x6; x7; x8; x9];
    
    vec = [[1; 0; 0.5], [1; 0; -0.5], [0; 1; 0.5], [0; 1; -0.5]]/2; % WLOG use eggbox topology
    [~,Pj,B,B_k] = parallelogramOrigami(2, vec);
    L = periodicityMatrix(B_k, 9, 3);

end

% Construct necessary tensors
if test < 5 % 2D case
    
    N = null(L_0); % null space matrix, size 2I by N where N = 2I-2(|B1|+|B2|+1)

    % construct homogeneous strain matrix X
    X = constructX(x);
    
    % compute initial moduli matrix C
    C_0 = effectiveC(x,X,N,B); % C is 3 by 3, P is 2I by 2I

else % 3D case (origami)

    N = null(L); % null space matrix for 3D unit cell, size 3I by N
    
    % construct homogeneous strain matrix X
    X = constructY(x, "translation");
    
    % compute initial moduli matrix C
    C_0 = membraneStiffness(x, Pj, L, 3, "translation");

end

% smallest eigenpair
[V,D] = eig(C_0);
[~,idx] = min(diag(D));
e_0 = V(:,idx); % strain eigenvector with smallest eigenvalue
s_0 = [1; 1; 0]; % initial homogeneous stretch tensor S = I in Voigt not.

% Initialization checks
fprintf("Initial C eigenvalues:\n");
disp(eig(C_0))

if test<5, K_0 = barStiffness(x,B); GH_0 = N'*K_0*N;
else, [~,~,GH_0] = membraneStiffness(x,Pj,L,3,"translation"); end
fprintf("Initial GH matrix eigenvalues:\n");
disp(eig(GH_0))

% if trace of strain tensor < 0, flip e_0 to enforce expansion
%if (e_0(1) + e_0(2)) < 0 || ~reverse
if reverse
    e_0 = -e_0;
end

state_0 = [e_0; s_0; x]; % e_0 = state_0(1:3), s = state(4:6), x = state_0(7:end)
% initial state variable checks
%{
e = e_0;
y = x;
[C,P] = effectiveC(y,X,N,B); % C is 3 by 3, P is 2I by 2I
D = eig(C);
lambda_min = min(D);
u = P*X*e; % column vector length 2I or 3I (origami)
[~, min_idx] = min(abs(e)); % component of e w/ smallest magnitude
e_aux = zeros(3,1);
e_aux(min_idx) = 1; % auxiliary vector along that axis
v1 = cross(e, e_aux);
v2 = cross(e, v1);
v1 = v1 / norm(v1);
v2 = v2/norm(v2);
P_0 = [v1, v2]; % 3 by 2 orthonormal matrix
C_0 = P_0'*(C - lambda_min*eye(3))*P_0; % 2 by 2 matrix
Kdot = zeros(length(y));
if test<5, chi = construct_chi2D(y);
else, chi = construct_chi(y);
end
for b = 1:size(pars.B,1)
    i = B(b,1);
    j = B(b,2);
    chi_ij = chi{i} - chi{j};
    y_ij = chi_ij*y;
    n = y_ij/norm(y_ij);
    Pij = eye(length(n)) - n*n'; % dim by dim projector matrix
    Kdot = Kdot + chi_ij'*(Pij*chi_ij*u/norm(y_ij)*n'+n*(Pij*chi_ij*u/norm(y_ij))')*chi_ij;
end
DuC = X'*P'*Kdot*P*X;
f = -P_0*(C_0\P_0')*DuC*e;
%}
% parameters for ODE function
pars.X = X;
pars.B = B;
pars.N = N;

pars.testNum = test;
if test == 1, pars.l_0 = [x(3:4), x(7:8)]; pars.boundary_pairs = [1 2; 1 4]; end
if test == 2, pars.l_0 = [x(5:6), x(13:14)]; pars.boundary_pairs = [1 3; 1 7]; end
if test == 3, pars.l_0 = [x(9:10), x(19:20)-x(3:4)]; pars.boundary_pairs = [1 5; 2 10]; end
if test == 4, pars.l_0 = [x(5:6), x(11:12)]; pars.boundary_pairs = [1 3; 1 6]; end
if test >= 5, pars.L = L; pars.Pj = Pj; pars.dim = 3; pars.group = "translation"; end

%% Integrate the system

%options = odeset('Events', @(t,state) combinedEvents(t, state, pars));
options = odeset('Events', @(t,state) GHEvent(t, state, pars), ...
    'RelTol',1e-7, ...
    'AbsTol',1e-9);
[t,state,te,ye,ie] = ode45(@(t,state) mechanismODE(t,state,pars),...
                           [0,T],state_0,options);
disp("Integration complete.")

%% Post processing

eigResidual = zeros(length(t),1);
rayleighResidual = zeros(length(t),1);
normResidual = zeros(length(t),1);
eigResidual_eig = zeros(length(t),1);
eVectorError = zeros(length(t),1);

eigs = zeros(3,1,length(t));
GH_eig = zeros(length(t),1);

% mechanism animation
figure
axis equal
if test>4, view(3)
end
hold on
frameDelay = 1/60; % 0.05 ~20 FPS
for k = 1:length(t)

    y = state(k,7:end)';
    
    
    if test>1 && test<5
    K = barStiffness(y,B);
    GH_matrix = N'*K*N; % must be invertible; otherwise GH-mode is present
    GH_matrix = (GH_matrix+GH_matrix')/2; % symmetricize
    GH_matrix(abs(GH_matrix)<tol) = 0; % round off numerics
    GH_eigs = eig(GH_matrix);
    GH_eig(k) = min(GH_eigs);
    end
    
    if test < 5
        C = effectiveC(y,X,N,B);
    else
        C = membraneStiffness(y, Pj, L, pars.dim, pars.group);
    end
    eigs(:,:,k) = eig(C);

    if test<5, Y = reshape(y,2,[])';
    else,      Y = reshape(y,3,[])'; end

    cla

    [V,D] = eig(C);
    %lambda_min = min(eig(C));
    e = state(k,1:3)';
    [lambda_min,idx] = min(diag(D));
    
    e_eig = V(:,idx);
    
    % Fix arbitrary sign of eig()
    if dot(e,e_eig) < 0
        e_eig = -e_eig;
    end

    eigResidual(k) = norm(C*e - lambda_min*e);
    rayleighResidual(k) = abs(e'*C*e - lambda_min);
    normResidual(k) = abs(norm(e)-1);
    eigResidual_eig(k) = norm(C*e_eig - lambda_min*e_eig);
    eVectorError(k) = norm(e-e_eig);
    if plotMechanism
    for b = 1:size(B,1)
        i = B(b,1);
        j = B(b,2);
        if test<5
            plot(Y([i j],1),Y([i j],2),'b-','LineWidth',2);
            if tessellate
                l1_1 = Y(pars.boundary_pairs(1,2),1)-Y(pars.boundary_pairs(1,1),1);
                l1_2 = Y(pars.boundary_pairs(1,2),2)-Y(pars.boundary_pairs(1,1),2);
                l2_1 = Y(pars.boundary_pairs(2,2),1)-Y(pars.boundary_pairs(2,1),1);
                l2_2 = Y(pars.boundary_pairs(2,2),2)-Y(pars.boundary_pairs(2,1),2);
                
                plot(Y([i j],1)+l1_1,Y([i j],2)+l1_2,'b-','LineWidth',2); % right tessellation
                plot(Y([i j],1)+l2_1,Y([i j],2)+l2_2,'b-','LineWidth',2); % top tessellation
                plot(Y([i j],1)+l1_1+l2_1,Y([i j],2)+l1_2+l2_2,'b-','LineWidth',2); % diagonal tessellation           
            end
        else, plot3(Y([i j],1),Y([i j],2), Y([i j],3),'b-','LineWidth',2); end
    end
    
    if test<5, plot(Y(:,1),Y(:,2),'ro','MarkerFaceColor','r');
    else, plot3(Y(:,1),Y(:,2), Y(:,3),'ro','MarkerFaceColor','r'); end

    for n = 1:size(Y,1)
        if test<5, text(Y(n,1) + 0.05, Y(n,2) + 0.03, num2str(n), ...
             'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
        else, text(Y(n,1) + 0.05, Y(n,2) + 0.03, Y(n,3) + 0.03, num2str(n), ...
             'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
        end
    end

    title(sprintf('t = %.3f',t(k)));
    drawnow;
    pause(frameDelay);
    end
end
hold off

% mechanism snapshots
if snapshotPlots
snapshot_coords = zeros(length(x), 5);
for num_snapshot=0:4
    if num_snapshot==0, k = 1; % initialization state
    elseif num_snapshot == 4, k = length(t); % final state
    else, k = num_snapshot*floor(length(t)/4); end % 3 intermediate states
    
    %if num_snapshot==3, k = k-floor(length(t)/8); end %

    y = state(k,7:end)';
    snapshot_coords(:,num_snapshot+1) = y;
    if test<5, Y = reshape(y,2,[])';
    else,      Y = reshape(y,3,[])'; end

    figure
    axis equal
    hold on
    if test>4, view(3)
    end

    for b = 1:size(B,1)
        i = B(b,1);
        j = B(b,2);
        if test<5, plot(Y([i j],1),Y([i j],2),'b-','LineWidth',2);
        else, plot3(Y([i j],1),Y([i j],2), Y([i j],3),'b-','LineWidth',2); end
    end
    
    if test<5, plot(Y(:,1),Y(:,2),'ro','MarkerFaceColor','r');
    else, plot3(Y(:,1),Y(:,2), Y(:,3),'ro','MarkerFaceColor','r'); end

    for n = 1:size(Y,1)
        if test<5, text(Y(n,1) + 0.03, Y(n,2) + 0.03, num2str(n), ...
             'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
        else, text(Y(n,1) + 0.05, Y(n,2) + 0.03, Y(n,3) + 0.03, num2str(n), ...
             'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
        end
    end
    title(sprintf('t = %.3f',t(k)));
    hold off
end
end

% eigenvalue plot
if eigPlot
zero_mode = squeeze(eigs(1,1,:)); % mechanism mode
soft_mode = squeeze(eigs(2,1,:));
stiff_mode = squeeze(eigs(3,1,:));
figure
hold on
plot(t, zero_mode, 'ko', 'MarkerFaceColor', 'k') 
plot(t, soft_mode, 'bo', 'MarkerFaceColor', 'b')
%plot(t, stiff_mode, 'ro', 'MarkerFaceColor', 'r')
xlabel("Time")
ylabel("Principal Membrane Stiffness")
%legend(["Zero mode", "Soft mode", "Stiff mode"])
hold off

figure
hold on
semilogy(t, GH_eig, 'ko', 'MarkerFaceColor', 'k')
%title('2x2 Square Lattice GH Matrix Min Eigenvalue')
xlabel("Time")
ylabel("\lambda_m (N^T KN)")
ylim([10^-3 1])
hold off
end

% energy check
if energyPlot || strainPlot

    E = zeros(length(t),1);
    strains = zeros(length(t), length(B));
    total_abs_strains = zeros(length(t));

    if test<5, X = reshape(x,2,[])';
    else X = reshape(x,3,[])';
    end

    for k=1:length(t)
        
        y = state(k,7:end)';
        if test<5, Y = reshape(y,2,[])';
        else Y = reshape(y,3,[])';
        end

        for b = 1:size(B,1)
            i = B(b,1);
            j = B(b,2);

            l0 = norm(X(i,:)-X(j,:));
            l = norm(Y(i,:)-Y(j,:));
            strains(k,b) = (l - l0)/l0;
            E(k) = E(k) + 0.5*(l - l0)^2;
        end
        total_abs_strains(k) = sum(abs(strains(k,:)));
    end
    if energyPlot
    figure
    semilogy(t, E)
    hold on
    xlabel("Time")
    ylabel("Energy")
    hold off
    end
    if strainPlot
    figure
    semilogy(t, total_abs_strains)
    hold on
    xlabel("Time")
    ylabel("Total absolute strain")
    hold off
    end
end

% residual check
if residualCheck
figure
semilogy(t,eigResidual,'LineWidth',2)
hold on
%semilogy(t, eigResidual_eig, 'LineWidth',2);
%semilogy(t, eVectorError, 'LineWidth',2);
%semilogy(t,rayleighResidual,'LineWidth',2)
semilogy(t,normResidual,'LineWidth',2)
%title('2x2 Square Lattice Error')
xlabel('Time')
ylabel('Residual')
ylim([10^-17 max(eigResidual)*10^2])
legend('||Ce-\lambda e||','||e||-1')
%  '||Ce_{eig}-\lambda e_{eig}||', 'e^TCe-\lambda_{min}', , '||e-e_{eig}||'
grid on
hold off
end

% homogeneous deformation tensor S check
%{
S_min_eig = zeros(length(t),1);
for k = 1:length(t)
    s = state(k,4:6)';
    S = [s(1)   s(3)/2;
         s(3)/2 s(2)];
    S_min_eig(k) = min(eig(S));
end
figure
hold on
plot(t,S_min_eig,'LineWidth',2)
xlabel('Time')
ylabel('\lambda_{min}(S)')
grid on
hold off
%}

[C0,P0] = effectiveC(x,X,N,B);

%% Functions

% ODE function, returns [edot; sdot; ydot] = [g; e; u]
function dq = mechanismODE(t,state,pars)

e = state(1:3);
s = state(4:6);
y = state(7:end);

e = e/norm(e);

% compute initial moduli matrix C
if pars.testNum < 5
    [C,P] = effectiveC(y,pars.X,pars.N,pars.B); % C is 3 by 3, P is 2I by 2I
else
    [C,P] = membraneStiffness(y, pars.Pj, pars.L, pars.dim, pars.group);
    %X_current = constructY(y,pars.group);
end
% smallest eigenvalue
D = eig(C);
lambda_min = min(D);

% compute displacement vector u = ydot
u = P*pars.X*e; % column vector length 2I or 3I (origami)
% store strain vector e
sdot = e;

% compute tangent space projector Ptilde
[~, min_idx] = min(abs(e)); % component of e w/ smallest magnitude
e_aux = zeros(3,1);
e_aux(min_idx) = 1; % auxiliary vector along that axis
v1 = cross(e, e_aux);
v2 = cross(e, v1);
v1 = v1 / norm(v1);
v2 = v2/norm(v2);
P_0 = [v1, v2]; % 3 by 2 orthonormal matrix

% compute edot = g = Pt*lam2:
% first the LHS for lam2
C_0 = P_0'*(C - lambda_min*eye(3))*P_0; % 2 by 2 matrix
if rcond(C_0) < 1e-12
    warning('Ill-conditioned tangent operator at t = %.6f\n', t)
end
% next the RHS for lam2, for which we first need DuC
% compute derivative of C using derivative of K
Kdot = zeros(length(y));% Eq. (80)
if pars.testNum<5, chi = construct_chi2D(y);
else, chi = construct_chi(y);
end
for b = 1:size(pars.B,1) % goes down the list of index pairs
    i = pars.B(b,1); % first index
    j = pars.B(b,2); % second index
    chi_ij = chi{i} - chi{j};
    y_ij = chi_ij*y;
    n = y_ij/norm(y_ij);
    Pij = eye(length(n)) - n*n'; % dim by dim projector matrix
    Kdot = Kdot + chi_ij'*(Pij*chi_ij*u/norm(y_ij)*n'+n*(Pij*chi_ij*u/norm(y_ij))')*chi_ij;
end
DuC = pars.X'*P'*Kdot*P*pars.X;
f = -P_0*(C_0\P_0')*DuC*e;

%f_scale = floor(log10(norm(f))); % order of magnitude of f
%gamma = 10^(f_scale); % damping factor
%damping = gamma*(C*e - (e'*C*e)*e);
%f = f-damping;

dq = [f; sdot; u];

end

function K = barStiffness(y, B)

    tol = 10^(-15);

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

function [C,P] = effectiveC(y,X,N,B)

    tol = 10^(-15);
    
    % compute stiffness matrix K
    K = barStiffness(y,B);
    
    % compute projection matrix P
    GH_matrix = N'*K*N; % must be invertible; otherwise GH-mode is present
    GH_matrix = (GH_matrix+GH_matrix')/2; % symmetricize
    P = eye(length(y)) - N*((GH_matrix)\(N'*K)); % 2I by 2I matrix
    %P = eye(length(y)) - N*pinv(GH_matrix)*N'*K; % 2I by 2I matrix
    
    % compute moduli matrix C
    C = X'*P'*K*P*X; % 3 by 3 symmetric positive semi-definite matrix
    C = (C+C')/2; % symmetricize
end

function [value,isterminal,direction] = GHEvent(t,state,pars)

    tol = 10^(-12);
    
    y = state(7:end);
    
    if pars.testNum<5
        K = barStiffness(y,pars.B);
        GH = pars.N'*K*pars.N;
        GH = (GH+GH')/2; % symmetricize
        GH(abs(GH)<tol) = 0; % round off numerics
    else, [~,~,GH] = membraneStiffness(y,pars.Pj,pars.L,pars.dim,pars.group);
    end
    
    eigs = eig(GH);
    min_eig = min(eigs);
    value = min_eig - tol;

    if value < 0
        num_GH = length(eigs(abs(eigs)<tol));
        fprintf('%d GH modes detected at t = %.6f\n', num_GH, t);
    end
    
    isterminal = 1;              % terminate integration
    direction = -1;              % only detect decreasing through threshold

end

function [value, isterminal, direction] = PSDEvent(t, state, pars)

    tol = 10^(-12);

    e = state(1:3);
    s = state(4:6);
    y = state(7:end);    

    if pars.testNum < 5

        i1prime = pars.boundary_pairs(1,2);
        i1 = pars.boundary_pairs(1,1);
        i2prime = pars.boundary_pairs(2,2);
        i2 = pars.boundary_pairs(2,1);

        l1 = y(2*i1prime-1:2*i1prime) - y(2*i1-1:2*i1);
        l2 = y(2*i2prime-1:2*i2prime) - y(2*i2-1:2*i2);
        S = [l1, l2] / pars.l_0;
        value = det(S) - tol;

    else
        
        S = [s(1)   s(3)/2;
             s(3)/2 s(2)];
        value = min(eig(S)) - tol;

    end
    
    if value <= 0
        fprintf('PSD threshold crossed at t = %.6f\n', t);
    end
    isterminal = 1;
    direction = -1;
end

function [value, isterminal, direction] = combinedEvents(t, state, pars)
    [val_GH, term_GH, dir_GH] = GHEvent(t, state, pars);
    [val_PSD, term_PSD, dir_PSD] = PSDEvent(t, state, pars);
    
    % concatenate results
    value      = [val_GH; val_PSD];
    isterminal = [term_GH; term_PSD];
    direction  = [dir_GH; dir_PSD];
end
