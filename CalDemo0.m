% Use fem calculate 2d integral
% To simplify the program, we do the integrals on the same domain

%% 0. Quadrature test
clear
U1.Degree = 1;
U1.N = 4;
U1.Nt = 1;
U1.Node = [0,1; 1,0; 0,0];
U1.TC = [1,2,3];
U1.Tri = [1,2,3];
K1 = symBilinear(U1, 'mass');
u = @(x, y) x;
v = @(x, y) y;

% Get Nodal info
U = u(U1.Node(:,1), U1.Node(:,2));
V = v(U1.Node(:,1), U1.Node(:,2));

fprintf("0\n")
fprintf("The FEM approximation value is %.6f\n", U'*K1*V);
fprintf("The true value is %.6f\n", integral2(@(x,y) x.*y,0,1,0,@(x)1-x));

U2.Degree = 2;
U2.N = 3;
U2.Nt = 1;
U2.Node = [0,1; 0,0; 1,0; 0, 1/2; 1/2, 0; 1/2, 1/2];
U2.TC = [1,2,3,4,5,6];
U2.Tri = [1,4,6;4,2,5;4,5,6;6,5,3];
K2 = symBilinear(U2, 'mass');
u = @(x, y) x.^2;
v = @(x, y) x.*y;


% Get Nodal info
U = u(U2.Node(:,1), U2.Node(:,2));
V = v(U2.Node(:,1), U2.Node(:,2));

fprintf("The FEM approximation value is %.6f\n", U'*K2*V);
fprintf("The true value is %.6f\n", integral2(@(x,y) x.^3.*y,0,1,0,@(x)1-x));


%% Intergration Part
% Define the domain and generate mesh
clear
N = 4;
Nx = N; Ny = N;
T = RecMesh([Nx, Ny], [1, 1], [0, 0]);

% Define a space for calculation
P = P1Fespace(T);

G = {@(x, y)0, @(x, y)0, @(x, y)0, @(x, y)0};
[BX, FreeB] = Freedomdefine(P, [1,1,1,1], G);
NB = sum(FreeB);
A = 0.2;
T.Node(FreeB, :) = T.Node(FreeB, :) + A*rand(NB, 2) - A/2;
ShowMesh(T);

P = P1Fespace(T);
U = P2Fespace(T);

% Generate the matrix
M1 = symBilinear(P, 'mass', []);
M2 = symBilinear(U, 'mass', []);
% Kux = symBilinear(U, 'dx');
Kux = nonsymBilinear(T, U, 'dx', U, 'dx', []);
Kuxpx = nonsymBilinear(T, U, 'dx', P, 'dx', []);
Kuxp = nonsymBilinear(T, U, 'dx', P, 'mass', []);


%% 1. Calculate (1, u)
u = @(x, y) x.^4;

% Get Nodal info
U1 = u(P.Node(:,1), P.Node(:,2));
V1 = ones(P.N,1);
U2 = u(U.Node(:,1), U.Node(:,2));
V2 = ones(U.N,1);

fprintf("1\n")
fprintf("The first order FEM approximation value is %.6f\n", V1'*M1*U1);
fprintf("The second order FEM approximation value is %.6f\n", V2'*M2*U2);
fprintf("The true value is %.6f\n", integral2(u,0,1,0,1));

%% 2. Calculate (u, v)
u = @(x, y) x;
v = @(x, y) y.^2;

% Get Nodal info
U1 = u(P.Node(:,1), P.Node(:,2));
V1 = v(P.Node(:,1), P.Node(:,2));
U2 = u(U.Node(:,1), U.Node(:,2));
V2 = v(U.Node(:,1), U.Node(:,2));

fprintf("2\n")
fprintf("The first order FEM approximation value is %.6f\n", V1'*M1*U1);
fprintf("The second order FEM approximation value is %.6f\n", V2'*M2*U2);
fprintf("The true value is %.6f\n", integral2(@(x,y) u(x,y).*v(x,y),0,1,0,1));

%% 3. Calculate (dx(u), dx(v))
u = @(x, y) x.*y;
v = @(x, y) x.^2;

% Get Nodal info
U1 = u(U.Node(:,1), U.Node(:,2));
V1 = v(U.Node(:,1), U.Node(:,2));
V2 = v(P.Node(:,1), P.Node(:,2));

fprintf("3\n")
fprintf("The FEM approximation value is %.6f\n", U1'*Kux*V1);
fprintf("The FEM approximation value is %.6f\n", U1'*Kuxpx*V2);
fprintf("The true value is %.6f\n", integral2(@(x,y) 2*x.*y,0,1,0,1));

%% 4. Calculate (dx(u), v)
u = @(x, y) x.^2;
v = @(x, y) y;

% Get Nodal info
U1 = u(U.Node(:,1), U.Node(:,2));
V1 = v(P.Node(:,1), P.Node(:,2));

fprintf("4\n")
fprintf("The FEM approximation value is %.6f\n", U1'*Kuxp*V1);
fprintf("The true value is %.6f\n", integral2(@(x,y) 2*x.*y,0,1,0,1));

%% 5.Test of numP2Bilinear
u = @(x, y) x.^2;
f = @(x, y) x.*y;
v = @(x, y) y;

fnk = f(U.Node(:,1), U.Node(:,2));
Mat1 = numP2Bilinear(U, fnk, "dx", "mass");
Mat2 = numP2Bilinear(U, fnk, "mass", "dx");
U1 = u(U.Node(:,1), U.Node(:,2));
V1 = v(U.Node(:,1), U.Node(:,2));
fprintf("5\n")
fprintf("The FEM approximation value is %.6f\n", V1'*Mat1*U1);
fprintf("The true value is %.6f\n", integral2(@(x,y) 2*x.^2.*y.^2,0,1,0,1));
fprintf("The FEM approximation value is %.6f\n", V1'*Mat2*U1);
fprintf("The true value is %.6f\n", integral2(@(x,y) x.^2.*y.^2,0,1,0,1));


%% 6. Test of Bilinear1d
u = @(x, y) x.^2.*y;
v = @(x, y) y;

U1 = u(P.Node(:, 1), P.Node(:, 2));
V1 = v(P.Node(:, 1), P.Node(:, 2));
U2 = u(U.Node(:, 1), U.Node(:, 2));
V2 = v(U.Node(:, 1), U.Node(:, 2));

M1 = Bilinear1d(T, 3, P, "mass", P, "mass");
M2 = Bilinear1d(T, 3, U, "mass", U, "mass");
Mdx1 = Bilinear1d(T, 3, P, "dx", P, "mass");
Mdx2 = Bilinear1d(T, 3, U, "dx", U, "mass");
M3 = Bilinear1d(T, 3, U, "mass", P, "mass");
fprintf("6\n");
fprintf("The FEM approximation value is %.6f\n", V1'*M1*U1);
fprintf("The FEM approximation value is %.6f\n", U1'*M1*V1);
fprintf("The FEM approximation value is %.6f\n", U2'*M2*V2);
fprintf("The FEM approximation value is %.6f\n", U2'*M3*V1);
fprintf("The true value is %.6f\n", integral(@(x) x.^2, 0,1));
fprintf("The FEM approximation value is %.6f\n", U1'*Mdx1*V1);
fprintf("The FEM approximation value is %.6f\n", U2'*Mdx2*V2);
fprintf("The true value is %.6f\n", integral(@(x) 2*x, 0,1));

