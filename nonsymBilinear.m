function K = nonsymBilinear(U1, type1, U2, type2, fnk)
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
if degree == 1
    P = [ones(size(Quads2d.px)), Quads2d.px, Quads2d.py];
    Px = [zeros(size(Quads2d.px)), ones(size(Quads2d.px)), zeros(size(Quads2d.px))];
    Py = [zeros(size(Quads2d.px)), zeros(size(Quads2d.px)), ones(size(Quads2d.px))];
    rcord = [0, 0; 1, 0; 0, 1];
    C = [ones(size(rcord,1),1), rcord];
end
if degree == 2
    P = [ones(size(Quads2d.px)), Quads2d.px, Quads2d.py,...
        Quads2d.px.^2, Quads2d.px.*Quads2d.py, Quads2d.py.^2];
    Px = [zeros(size(Quads2d.px)), ones(size(Quads2d.px)), zeros(size(Quads2d.px)),...
        2*Quads2d.px, Quads2d.py, zeros(size(Quads2d.px))];
    Py = [zeros(size(Quads2d.px)), zeros(size(Quads2d.px)), ones(size(Quads2d.px)),...
        zeros(size(Quads2d.px)), Quads2d.px, 2*Quads2d.py];
    rcord = [0, 0; 1, 0; 0, 1; 1/2, 0; 1/2, 1/2; 0, 1/2];
    C = [ones(size(rcord,1),1), rcord, rcord(:,1).^2, rcord(:,1).*rcord(:,2), rcord(:,2).^2];
end
% M = inv(C); M is the coefficients matrix of basis funtions
if type1 == "mass"
    Ical1 = P/C;
else
    Ixcal1 = Px/C;
    Iycal1 = Py/C;
end
if type2 == "mass"
    Ical2 = P/C;
else
    Ixcal2 = Px/C;
    Iycal2 = Py/C;
end
indk = 1;
n1 = size(U1.TC(1, :), 2);
n2 = size(U2.TC(1, :), 2);
tmp = ceil(U1.Nt*n1*n2);
indi = ones(tmp ,1);
indj = ones(tmp, 1);
valuek = zeros(tmp, 1);
for i = 1:U1.Nt
    ind = U1.TC(i, :);
    cord = U1.Node(ind, :);
    A = (cord(2:3, :) - repmat(cord(1, :),2,1))';
    Area = abs(det(A));
    J = inv(A);
    for s = 1:n1
        for r = 1:n2
            indi(indk, 1) = ind(s);
            indj(indk, 1) = ind(r);
            switch type1
                case {"mass"}; I1 = Ical1(:, s);
                case {"dx"}; I1 = ([Ixcal1(:, s), Iycal1(:, s)]*J(:, 1));
                case {"dy"}; I1 = ([Ixcal1(:, s), Iycal1(:, s)]*J(:, 2));
                case {"nabla"}; I1 = ([Ixcal1(:, s), Iycal1(:, s)]/A);
            end
            switch type2
                case {"mass"}; I2 = Ical2(:, s);
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