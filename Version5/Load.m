function F = Load(U, f)
% This function assembles the Load Vector without Neumann boundary
% Condition

if nargin < 2; f = @(x, y) 1; end
degree = U.Degree;
Quads2d = LoadQuad2d(2*degree);
% Tr = RefInfo(U.Degree, Quads2d);
% Ical = Tr.P/Tr.C;
if degree == 1
    P = [ones(size(Quads2d.px)), Quads2d.px, Quads2d.py];
    rcord = [0, 0; 1, 0; 0, 1];
    C = [ones(size(rcord,1), 1), rcord];
end
if degree == 2
    P = [ones(size(Quads2d.px)), Quads2d.px, Quads2d.py,...
        Quads2d.px.^2, Quads2d.px.*Quads2d.py, Quads2d.py.^2];
    rcord = [0, 0; 1, 0; 0, 1; 1/2, 0; 1/2, 1/2; 0, 1/2];
    C = [ones(size(rcord,1), 1), rcord, rcord(:,1).^2, rcord(:,1).*rcord(:,2), rcord(:,2).^2];
end
Ical = P/C;
n = size(U.TC(1, :), 2);
N = U.Nt*n;
indi = ones(N, 1);
indj = ones(N, 1);
valueF = zeros(N, 1);
indk = 1;

for i = 1:U.Nt
    ind = U.TC(i, :);
    cord = U.Node(ind, :);
    A = (cord(2:3, :) - repmat(cord(1, :),2,1))';
    Area = abs(det(A));
    for s = 1:n
        cordf = repmat(cord(1, :), size(Quads2d.px, 1), 1) + [Quads2d.px, Quads2d.py]*A';
        indi(indk) = ind(s);
        valueF(indk) = Quads2d.w*(f(cordf(:,1), cordf(:,2)).*Ical(:, s))/2*Area;
        indk = indk + 1;
    end
end
F = sparse(indi, indj, valueF);
end