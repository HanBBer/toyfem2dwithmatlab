%G = {@(x,y)0, @(x,y)0, @(x,y)0, @(x,y)0};
u = @(x, y) exp(sin(x)+cos(y));
f = @(x, y) (cos(y)+sin(x)-cos(x).^2+cos(y).^2).*exp(sin(x)+cos(y));
G = {u, u, u, u};

N = 2;


for i = 1:3
    N = N*2;
    nx = N; ny = N;
    T = RecMesh([nx, ny], [1, 1], [0, 0]);
    U = P1Fespace(T);
    K = symBilinear(U, 'nabla', []) + symBilinear(U, 'mass', []);
    [X, FNodeptr] = Freedomdefine(U, [1,1,1,1], G);
    F = Load(U, f);
    KK = K(FNodeptr, FNodeptr);
    FF = F(FNodeptr) - K(FNodeptr, ~FNodeptr)*X(~FNodeptr);
    X(FNodeptr) = KK\FF;
    % H1 norm
    err = sqrt( (X - u(U.Node(:,1),U.Node(:,2)))'*K*(X - u(U.Node(:,1),U.Node(:,2))) );
    % L2 norm
    %err = symBilinear(U, 'nabla', [])*(X - u(U.Node(:, 1),U.Node(:, 2)));
    fprintf('Error: \t%.4e \n',err)
    %trisurf(U.Tri, U.Node(:, 1), U.Node(:, 2), X - u(U.Node(:, 1),U.Node(:, 2)));
    if i>1
        fprintf('Rate: \t%.4f \n',log(errold/err)/log(2))
    end
    errold = err;
end