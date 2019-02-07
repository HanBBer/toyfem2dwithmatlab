function K = symBilinear(U, type, fnk)
% This function calculates the inner product of finite element function
% We implements this process element by element
% Here we do not treat any kind of boundary condition
% To get a sparse result efficiently, 
% 1. We introduce three vectors: indi indj valuek to record the sparse info
% 2. We do calculation on refference trigular
% 3. Since the matrix is symmetry, we do not calculate all the element,
%    but make a duplication to raise effciency
if nargin < 3; fnk = []; end
degree = U.Degree;
Quads2d = LoadQuad2d(2*degree);
Tr = RefInfo(U.Degree, Quads2d);
% if degree == 1
%     P = [ones(size(Quads2d.px)), Quads2d.px, Quads2d.py];
%     Px = [zeros(size(Quads2d.px)), ones(size(Quads2d.px)), zeros(size(Quads2d.px))];
%     Py = [zeros(size(Quads2d.px)), zeros(size(Quads2d.px)), ones(size(Quads2d.px))];
%     rcord = [0, 0; 1, 0; 0, 1];
%     C = [ones(size(rcord,1),1), rcord];
% end
% if degree == 2
%     P = [ones(size(Quads2d.px)), Quads2d.px, Quads2d.py,...
%         Quads2d.px.^2, Quads2d.px.*Quads2d.py, Quads2d.py.^2];
%     Px = [zeros(size(Quads2d.px)), ones(size(Quads2d.px)), zeros(size(Quads2d.px)),...
%         2*Quads2d.px, Quads2d.py, zeros(size(Quads2d.px))];
%     Py = [zeros(size(Quads2d.px)), zeros(size(Quads2d.px)), ones(size(Quads2d.px)),...
%         zeros(size(Quads2d.px)), Quads2d.px, 2*Quads2d.py];
%     rcord = [0, 0; 1, 0; 0, 1; 1/2, 0; 1/2, 1/2; 0, 1/2];
%     C = [ones(size(rcord,1),1), rcord, rcord(:,1).^2, rcord(:,1).*rcord(:,2), rcord(:,2).^2];
% end
% M = inv(C); M is the coefficients matrix of basis funtions
if type == "mass"
    Ical = Tr.P/Tr.C;
else
    Ixcal = Tr.Px/Tr.C;
    Iycal = Tr.Py/Tr.C;
end
indk = 1;
n = size(U.TC(1, :), 2);
tmp = ceil(U.Nt*n*(n+1)/2);
indi = ones(tmp ,1);
indj = ones(tmp, 1);
valuek = zeros(tmp, 1);
for i = 1:U.Nt
    ind = U.TC(i, :);
    cord = U.Node(ind, :);
    A = (cord(2:3, :) - repmat(cord(1, :), 2,1))';
    Area = abs(det(A));
    J = inv(A);
    for s = 1:n
        for r = 1:s
            indi(indk, 1) = ind(s);
            indj(indk, 1) = ind(r);
            switch type
                case {"mass"}; I = Ical(  :, s).* Ical(:, r);
                case {"dx"}; I = ([Ixcal(:, s), Iycal(:, s)]*J(:, 1)).*([Ixcal(:, r), Iycal(:, r)]*J(:, 1));
                case {"dy"}; I = ([Ixcal(:, s), Iycal(:, s)]*J(:, 2)).*([Ixcal(:, r), Iycal(:, r)]*J(:, 2));
                case {"nabla"}; I = ([Ixcal(:, s), Iycal(:, s)]/A).*([Ixcal(:, r), Iycal(:, r)]/A);
            end
            if ~isempty(fnk)
                cordf = repmat(cord(1, :), n,1) + [Quads2d.px, Quads2d.py]*A';
                valuek(indk, 1) = sum(Quads2d.w*(fnk(cordf(:,1), cordf(:,2)).*I))/2*Area;
            else
                valuek(indk, 1) = sum(Quads2d.w*I)/2*Area;
            end
            indk = indk + 1;
        end
    end
end
K = sparse(indi, indj, valuek);
K = K + K' - diag(diag(K));
end