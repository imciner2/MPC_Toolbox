%% Create the physical system for the test
A = [0.7, -0.1, 0.0, 0.0;
     0.2, -0.5, 0.1, 0.0;
     0.0,  0.1, 0.1, 0.0;
     0.5,  0.0, 0.5, 0.5];
 
B = [0.0, 0.1;
     0.1, 1.0;
     0.1, 0.0;
     0.0, 0.0];

[n, m] = size(B);
 
% Create the Q matrix
Q = eye(n);
Q(1,1) = 10;
Q(2,2) = 20;
Q(3,3) = 30;
Q(4,4) = 40;

% Create the R matrix
R = [10,  0;
      0, 20];

P = dlyap(A', Q);


%% Create the constraints
ulb = [1;
      -1];

E  = kron(eye(m), ulb );
Ce = kron( ones(m,1), [ 0.5; 0.5] );

D  = kron(eye(n), ulb );
Cd = kron( ones(n,1), [ 0.5; 0.5] );

%% Create some constants
% Prediction horizon
N = 20;

% Initial condition
xStart = [0.1;
          0.1;
          0.1;
          0.1];

% Error terms
eps_g  = 1e-4;    % Constraint satisfaction error

% How closely the solutions should match
TOL = 1e-4;

%% Generate the sparse problem matrices
Hs = sparse.gen_cost(N, Q, R, P);
Aeq = sparse.gen_con_eq(N, A, B);
[Aineq, bineq] = sparse.gen_con_ineq(N, n, m, E, Ce, D, Cd);


%% Generate the condensed problem matrices
[Hp, Jp]  = condensed.primal.gen_cost(N, A, B, Q, R, P);
[G, F, g] = condensed.primal.gen_con_ineq(N, A, B, E, Ce, D, Cd);
[Hd, Jd]  = condensed.dual.gen_cost( Hp, Jp, G, F );

q  = Jp*xStart;
qd = Jd*xStart;


%% Generate options for quadprog to hide output
opts = optimoptions('quadprog');
opts.Display = 'off';


%% Test the sparse primal formulation using quadprog
beq = [xStart;
       zeros(n*N, 1)];
zHs = zeros( length(Hs), 1);

z = quadprog(Hs, zHs, Aineq, bineq, Aeq, beq, [], [], [], opts );

% Extract the control inputs
u_sparse_quad = [];
for (i=0:1:N-1)
    sI = (i+1)*n + i*m + 1;
    eI = (i+1)*n + (i+1)*m;
    u_sparse_quad = [u_sparse_quad;
                     z(sI:eI)];
end


%% Test the condensed primal formulation using quadprog
u_cond_quad = quadprog(Hp, q, G, (F*xStart + g), [], [], [], [], [], opts);


% Make sure the values are the same
e = u_cond_quad - u_sparse_quad;
if ( e > TOL )
    warning('Condensed primal formulation quadprog result is out of tolerance');
else
    disp('Condensed primal formulation quadprog result is inside the tolerance');
end


%% Test the dual formulation using quadprog
% Run the algorithm
y_quad = quadprog(Hd, (Jd*xStart + g), [], [], [], [], zeros(240, 1), inf(240, 1), [], opts);

% Recover the primal solution
u_y_quad = condensed.dual.recover_primal(y_quad, xStart, Hp, Jp, G);

% Make sure the values are the same
e = u_cond_quad - u_y_quad;
if ( e > TOL )
    warning('Condensed dual formulation quadprog result is out of tolerance');
else
    disp('Condensed dual formulation quadprog result is inside the tolerance');
end


%% Test the DGP algorithm
% Run the algorithm
[y_dgp, iters_dgp] = dgp_quad( Hp, q, G, g, eigs(Hd, 1, 'LM'), eps_g );

% Recover the primal solution
u_dgp = condensed.dual.recover_primal(y_dgp, xStart, Hp, Jp, G);

% Make sure the values are the same
e = u_cond_quad - u_dgp;
if ( e > TOL )
    warning('DGP result is out of tolerance.');
else
    disp('DGP result is inside the tolerance');
end


%% Test the FGM algorithm
% The projection operation
projOp = @(x) min( Ce(1), max(-Ce(2), x) );

% The suboptimality level
delta = 0.001;
subOpt_fgm = min(eig(Hp))/2 * delta^2/(norm(B,2)^2);

% Run the algorithm
[u_fgm, iters_fgm] = fgm_quad( Hp, q, projOp, eigs(Hp, 1, 'LM'), eigs(Hp, 1, 'SM'), zeros(m*N), 'Gradient', subOpt_fgm, 'Constant' );

% Make sure the values are the same
e = u_cond_quad - u_fgm;
if ( e > TOL )
    warning('FGM result is out of tolerance.');
else
    disp('FGM result is inside the tolerance');
end