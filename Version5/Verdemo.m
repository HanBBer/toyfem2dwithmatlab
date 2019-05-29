%dbstop if error;
debugflag = [0, 0, 0, 0, 0, 1];
%            1  2  3  4  5  6

% This is a script for testing

N = 1;
Nx = N; Ny = N;
T = RecMesh([Nx, Ny], [1,1], [0,0]);
if debugflag(1)
    figure(1); ShowMesh(T, [1,1,1]);
end

T2 = Refine(T);
if debugflag(2)
    figure(2); ShowMesh(T2, [1,1,1]);
end

% Verify the validity of Refine.m
T3 = Refine(T2);
if debugflag(3)
    figure(3); ShowMesh(T2, [1,1,1]);
end

% Test the P2 Space Generator
U2 = P2Fespace(T2);
if debugflag(4)
    figure(4); triplot(U2.Tri, U2.Node(:,1), U2.Node(:,2));
    box off; set(gca, 'XTick', [], 'YTick', []);
end

% Test LoadQuad2d
if debugflag(5)
    fprintf("Function \t Integral \t Quads \n")
    for i = 1:5
        j = randi([0,i]);
        f = @(x, y) x.^j.*y.^(i-j);
        Quads = LoadQuad2d(i);
        I1 = Quads.w2d*f(Quads.px2d, Quads.py2d)/2;
        I2 = integral2(f,0,1,0,@(x)1-x);
        fprintf("x^%d*y^%d \t %.6f \t %.6f\n",...
        i, j, I2, I1);
    end
end

K = symBilinear(U2, 'nabla', []);
%K = symBilinear(U2, 'dx', []) + symBilinear(U2, 'dy', []);
%K = nonsymBilinear(T, U2, 'nabla', U2, 'nabla', []);
g0 = @(x,y) 0;
G = {g0, g0, g0, g0};
[X, FNodeptr] = Freedomdefine(U2, [1,1,1,1], G);
F = Load(U2); 

KK = K(FNodeptr, FNodeptr);
FF = F(FNodeptr) - K(FNodeptr, ~FNodeptr)*X(~FNodeptr);
X(FNodeptr) = KK\FF;
if debugflag(6)
    figure;
    trisurf(U2.Tri, U2.Node(:, 1), U2.Node(:, 2), X);
end
