# Dynamical-system-for-constructing-cell-based-mechanisms-in-origami-and-kirigami-metamaterials
A cell-based mechanism in a periodic metamaterial is a smooth family of rigid deformations composed of a shape-changing homogeneous part along with a periodic correction that preserves the rigidity of the polygonal panels or bars in the setup. We implement a program to find these mechanisms that solve a dynamical system of the form

```math
\begin{bmatrix}
\dot{\mathbf{e}}(t) \\
\dot{\mathbf{y}}(t)
\end{bmatrix}
=
\begin{bmatrix}
\mathbf{f}(\mathbf{e}(t),\mathbf{y}(t)) \\
\mathbf{u}(\mathbf{e}(t),\mathbf{y}(t))
\end{bmatrix}
\quad \text{subject to} \quad
\begin{bmatrix}
\mathbf{e}(0) \\
\mathbf{y}(0)
\end{bmatrix}
=
\begin{bmatrix}
\mathbf{e}_0 \\
\mathbf{x}
\end{bmatrix}
```

where $\mathbf{y}(t)$ is an array listing the current positions of the vertices in a unit cell and $\mathbf{e}(t)$ is an array of the three in-plane components of the infinitesimal strain associated to an effective linear elastic
cell energy $\mathbf{e}\cdot\mathbf{C}(\mathbf{y})\mathbf{e}$. In particular, if there is exactly one shape-changing infinitesimal zero mode at the initial configuration $\mathbf{x}$, identified by the strain array $\mathbf{e}_0 \in \mathbb{S}^2$ such that $\mathbf{C}(\mathbf{x})\mathbf{e}_0 = \mathbf{0}$, and if there are no Guest-Hutchinson modes at $\mathbf{x}$, then the solution to this dynamical system constructs the mechanism motion whenever such a motion exists.

To run the program:
1. Initialize the following:
   - Nodal positions $\mathbf{x}\in\mathbb{R}^{2I}$ for rigid bar structures, or $\mathbf{x}\in\mathbb{R}^{3I}$ for rigid panel structures where $I$ is the number of nodes.
   - Bar connectivity $B\in\mathbb{Z}^{b\times 2}$ where $b$ is the number of bars or panel edges.
   - Periodicity pairs $B_k\in\mathbb{Z}^{(|B_1|+|B_2|)\times 2}$ where $B_k$ for $k=1,2$ is the set of index pairs related by periodic displacements $\mathbf{d}_k$.
   - Rigidity constraint matrix $L_0\in\mathbb{R}^{2(|B_1|+|B_2|+1)\times 2I}$ or $L\in\mathbb{R}^{3(|B_1|+|B_2|+1)\times 3I}$ constructed using function periodicityMatrix(B_k, num_nodes, dim).
   - (Rigid panel origami only): Panel labeling cell array $P_j$ where the $j$-th cell holds node indices of the $j$-th panel.
   - Parameters: integration time $T$ and boolean "reverse" to initialize the opposite mechanism motion.
2. Construct:
   - Null space matrix $N$ as
     ```matlab
     N = null(L_0);
     ```
     for 2D structures or
     ```matlab
     N = null(L);
     ```
     for 3D structures.
   - Mapping matrix $X$ or $Y$ using
     ```matlab
     X = constructX(x);
     ```
     for 2D structures or
     ```matlab
     Y = constructY(x,"group");
     ```
     for 3D structures.
   - Initial stiffness tensor $C_0$ using
     ```matlab
     C_0 = effectiveC(x,X,N,B);
     ```
     for rigid bar structures or
     ```matlab
     C_0 = membraneStiffness(x,Pj,L,dim,"group")
     ```
     for rigid panel structures.
3. Check initialization:
   - Confirm there are no Guest-Hutchinson modes before integrating by calling
     ```matlab
     K = barStiffness(x,B);
     GH_0 = N'*K*N;
     eig(GH_0);
     ```
     for rigid bar structures or
     ```matlab
     [~,~,GH0] = membraneStiffness(x,Pj,L,dim,"group")
     ```
     for rigid panel structures. If there is a zero eigenvalue you must choose a different initialization.
4. Prepare to integrate:
   - Use
     ```matlab
     state_0 = [e_0; x]
     ```
     and pass the parameters
     ```matlab
     pars.B = B; pars.N = N
     ```
     and for rigid bar structures:
     ```matlab
     pars.X = X, pars.l_0 = [l_1, l_2]
     ```
     for lattice vectors $\mathbf{l}_k$, and
     ```matlab
     pars.boundary_pairs = [i_1, i_1'; i_2, i_2']
     ```
     which indicate the corresponding node indices to the lattice vectors defined for pars.l_0. For rigid panel structures, instead pass:
     ```matlab
     pars.X = Y; pars.L = L; pars.Pj = Pj; pars.dim = dim; pars.group = "group"
     ```
5. Integrate by calling:
   ```matlab
   options = odeset('Events', @(t,state) combinedEvents(t, state, pars));
   [t,state,te,ye,ie] = ode45(@(t,state) mechanismODE(t,state,pars),[0,T],state_0,options);
   ```
