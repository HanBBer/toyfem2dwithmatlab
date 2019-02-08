function F = Load2(U, f, h)
% This function assembles the Load Vector
% It returns (f, v) - (h, v)_p
if nargin < 2; f = @(x, y) 1; end
hflag = [0,0,0,0];
if nargin == 3
    for i = 1:4
        hflag(i) = ~isempty(h{i});
    end
end
degree = U.Degree;
Quads2d = LoadQuad2d(2*degree);
Tr = RefInfo(U.Degree, Quads2d);
Ical = Tr.P/Tr.C;
N = 0;
if sum(hflag)
    Quads1d = LoadQuad1d(2*degree);
    nquads1d = size(Quads1d.px, 1);
    if degree == 1; Ical1d = [ones(nquads1d,1), Quads1d.px]/[1, -1; 1, 1]; end
    if degree == 2; Ical1d = [ones(nquads1d,1), Quads1d.px, (Quads1d.px).^2]/[1, -1, 1; 1, 0, 0; 1, 1, 1]; end
    % We need to estimate pre-allocated space for F
    if degree == 1; N = U.N+U.Ne-2*U.Nt; end
    if degree == 2; N = (U.N-U.Nt)/2*3; end
end
n = size(U.TC(1, :), 2);
N = U.Nt*n+N*sum(hflag);
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
    % Since the data structure is not so good to make the Neumann Boundary
    % easy-handled, we make an effort to realize the function
    % Note that the Node recorded in T.Tri is anticlockwise
    % So P1 case is fine
    % In P2 case U.TC has the order [T.Tri(i, 1:3), T.N+T.TrEg(i, 1:3)]
    % And T.TrEg(i, 1) is T.Tri(i, [1:2])
    % T.TrEg(i, 2) is T.Tri(i ,[2,3])  T.TrEg(i,3) is T.Tri(i. [3,1])
    % Which means the j-th edge is [j, j+3, mod(j+1,3)]
    % and the order is anticlockwise
    % so we can implemented NBC basiclly with the above-mentioned
    % consideration
    if sum(hflag) && sum(U.EgFlag(U.TrEg(i,:)))
        for j = 1:3
            egflag = U.EgFlag(U.TrEg(i, j));
            if egflag && hflag(egflag)
                g = h{egflag};
                if degree == 1
                    switch j
                        case{1}; ind = ind([1,2]);
                        case{2}; ind = ind([2,3]);
                        case{3}; ind = ind([3,1]);
                    end
                end
                if degree == 2
                    switch j
                        case{1}; ind = ind([1,4,2]);
                        case{2}; ind = ind([2,5,3]);
                        case{3}; ind = ind([3,6,1]);
                    end
                end
                egnode1 = U.Node(ind(1),:);
                egnode2 = U.Node(ind(end),:);
                len = norm(egnode1-egnode2,2);
                cordh = 1/2*repmat(egnode1+egnode2, nquads1d,1)+1/2*(egnode2-egnode1).*Quads1d.px;
                I = len/2*Ical1d'*(Quads1d.w'.*g(cordh(:, 1), cordh(:, 2)));
                for k = 1:length(ind)
                    indi(indk) = ind(k);
                    valueF(indk) = I(k);
                    indk = indk + 1;
                end
            end
        end
    end
end
F = sparse(indi, indj, valueF);
end