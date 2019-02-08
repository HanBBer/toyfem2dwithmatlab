function K = nonsymBilinear(T, U1, type1, U2, type2, fnk)
% This function calculates the inner product of finite element function
% We implements this process element by element
% Here we do not treat any kind of boundary condition
% To get a sparse result efficiently, 
% 1. We introduce three vectors: indi indj valuek to record the sparse info
% 2. We do calculation on refference trigular
% 3. Since the matrix is symmetry, we do not calculate all the element,
%    but make a duplication to raise effciency
if nargin < 5; fnk = []; end
degree = max(U1.Degree, U2.Degree);
Quads2d = LoadQuad2d(2*degree);
Tr1 = RefInfo(U1.Degree, Quads2d);
Tr2 = RefInfo(U2.Degree, Quads2d);
% M = inv(C); M is the coefficients matrix of basis funtions
if type1 == "mass"
    Ical1 = Tr1.P/Tr1.C;
else
    Ixcal1 = Tr1.Px/Tr1.C;
    Iycal1 = Tr1.Py/Tr1.C;
end
if type2 == "mass"
    Ical2 = Tr2.P/Tr2.C;
else
    Ixcal2 = Tr2.Px/Tr2.C;
    Iycal2 = Tr2.Py/Tr2.C;
end
indk = 1;
n1 = size(U1.TC(1, :), 2);
n2 = size(U2.TC(1, :), 2);
tmp = ceil(U1.Nt*n1*n2);
indi = ones(tmp ,1);
indj = ones(tmp, 1);
valuek = zeros(tmp, 1);
for i = 1:T.Nt
    ind = T.Tri(i, :);
    cord = T.Node(ind, :);
    A = (cord(2:3, :) - repmat(cord(1, :),2,1))';
    Area = abs(det(A));
    J = inv(A);
    ind1 = U1.TC(i, :);
    ind2 = U2.TC(i, :);
    for s = 1:n1
        for r = 1:n2
            indi(indk, 1) = ind1(s);
            indj(indk, 1) = ind2(r);
            switch type1
                case {"mass"}; I1 = Ical1(:, s);
                case {"dx"}; I1 = ([Ixcal1(:, s), Iycal1(:, s)]*J(:, 1));
                case {"dy"}; I1 = ([Ixcal1(:, s), Iycal1(:, s)]*J(:, 2));
                case {"nabla"}; I1 = ([Ixcal1(:, s), Iycal1(:, s)]/A);
            end
            switch type2
                case {"mass"}; I2 = Ical2(:, r);
                case {"dx"}; I2 = ([Ixcal2(:, r), Iycal2(:, r)]*J(:, 1));
                case {"dy"}; I2 = ([Ixcal2(:, r), Iycal2(:, r)]*J(:, 2));
                case {"nabla"}; I2 = ([Ixcal2(:, r), Iycal2(:, r)]/A);
            end
            if ~isempty(fnk)
                cordf = repmat(cord(1, :), n1,1) + [Quads2d.px, Quads2d.py]*A';
                valuek(indk, 1) = sum(Quads2d.w*(fnk(cordf(:,1), cordf(:,2)).*I1.*I2))/2*Area;
            else
                valuek(indk, 1) = sum(Quads2d.w*(I1.*I2))/2*Area;
            end
            indk = indk + 1;
        end
    end
end
K = sparse(indi, indj, valuek);

end