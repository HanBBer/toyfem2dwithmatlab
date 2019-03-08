function K = numP2Bilinear(U, fnk, typeU, typek)
% This function is designed intentionally for the useage in FSI problem
% it calculates the following bilinear matrax
% (fnk*dx(u), v) (dx(fnk)*u, v)
% where fnk/u/v are P2 FEspace

n = size(U.TC(1, :), 2);
tmp = U.Nt*n*n;
indi = ones(tmp ,1);
indj = ones(tmp, 1);
valuek = zeros(tmp, 1);
indk = 1;

Quads2d = LoadQuad2d(4);
Tr = RefInfo(2, Quads2d);
Ical = Tr.P/Tr.C;
Ixcal = Tr.Px/Tr.C;
Iycal = Tr.Py/Tr.C;


for i = 1:U.Nt
    ind = U.TC(i, 1:3);
    cord = U.Node(ind, :);
    A = (cord(2:3, :) - repmat(cord(1, :),2,1))';
    Area = abs(det(A));
    J = inv(A);
    ind = U.TC(i, :);
    
    switch typek
        case {"mass"}; I2 = Ical*fnk(ind);
        case {"dx"}; I2 = reshape( ([Ixcal(:), Iycal(:)]*J(:, 1)), 6, 6)*fnk(ind);
        case {"dy"}; I2 = reshape( ([Ixcal(:), Iycal(:)]*J(:, 2)), 6, 6)*fnk(ind);
    end
    
    for s = 1:n
        for r = 1:n
            indi(indk, 1) = ind(r);
            indj(indk, 1) = ind(s);
            switch typeU
                case {"mass"}; I1 = Ical(:, s);
                case {"dx"}; I1 = ([Ixcal(:, s), Iycal(:, s)]*J(:, 1));
                case {"dy"}; I1 = ([Ixcal(:, s), Iycal(:, s)]*J(:, 2));
            end
            
            I3 = Ical(:, r);
            valuek(indk, 1) = Quads2d.w*(I1.*I2.*I3)/2*Area;
            indk = indk + 1;
        end
    end
end
K = sparse(indi, indj, valuek);

end