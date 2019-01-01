function F = FemLinearLoad(T, Fd, f, h)
U = Fd.Space;
Nf = size(Fd.FNodePtrs, 1);
F = zeros(Nf, 1);
ns = size(T.(U).TC, 2);
[w, P, ~, ~, C01, C02] = LoadQuad();
Pcord = P(:,2:3);
if (T.(U1).Property == "P1");   C = C01;
else;  C = C02; end
for i = 1:T.Nt
    if ~isempty(f) || ~isempty(h)
        j = T.Tri(i, :);
        cord = T.Nodes(j,:);
        A = (cord(2:3,:) - repmat(cord(1,:),2,1))';
        Pcal = Pcord*A' + cord(1,:);
        Area = abs(det(A));
        jp = Fd.NodePtrs( T.(U).TC(i, :) );
    end
    
    if ~isempty(f)
        for s = 1:ns
            if jp(s) > 0
                F(s) = F(s) + w*f(Pcal(:,1), Pcal(:,2)).*(P*C(:, s))*Area;
            end
        end
    end
    
    if ~isempty(h)
        for s = 1:ns
            if jp(s) > 0
                F(s) = F(s);
            end
        end
    end
end
end