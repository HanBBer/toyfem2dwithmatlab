Lame = lambda*[Sx, Sxy; Sxy', Sy] + muS*[2*Sx+Sy, Sxy'; Sxy, 2*Sy+Sx];
FSIS = [[rhoS/dt*MS, sparse(NSu, NSu); sparse(NSu, NSu), rhoS/dt*MS], Lame;
    -[MS, sparse(NSu, NSu); sparse(NSu, NSu), MS], [MS/dt, sparse(NSu, NSu); sparse(NSu, NSu), MS/dt]];

FSISold = [[rhoS/dt*MS, sparse(NSu, NSu); sparse(NSu, NSu), rhoS/dt*MS], sparse(2*NSu, 2*NSu);
    sparse(2*NSu, 2*NSu), [MS/dt, sparse(NSu, NSu); sparse(NSu, NSu), MS/dt]];

G = {[], @(x, y) 0, [], @(x, y) 0};
[X, FNodeptr] = Freedomdefine(US, [0,1,0,1], G);

Free = [FNodeptr;FNodeptr;FNodeptr;FNodeptr];

Xold = XFSIold(2*NFu+NFp+1:end);

BigF = sparse(4*NSu, 1);
KK = FSIS(Free, Free);
XX = [X;X;X;X];

for t = 0:dt:2000*dt
    FF = BigF + FSISold*Xold;
    FF = FF(Free) - FSIS(Free, ~Free)*XX(~Free);
    XX(Free) = KK\FF;
    
    triplot(TS.Tri, TS.Node(:, 1)+XX(2*NSu+1:2*NSu+NSp), TS.Node(:, 2)+XX(3*NSu+1:3*NSu+NSp));
    axis equal
    pause(0.01)
    Xold = XX;
end