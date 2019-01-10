function F = FemLinearLoad(T, Fd, f, h)
U = Fd.Space;
Nf = size(Fd.FNodePtrs, 1);
F = zeros(Nf, 1);
ns = size(T.(U).TC, 2);
[w, P, ~, ~, C01, C02] = LoadQuad();
Pcord = P(:,2:3);
if (T.(U).Property == "P1");   C = C01;
elseif (T.(U).Property == "P2");  C = C02;
else; fprintf("error");
end
V = P*C;

for i = 1:T.Nt
    if ~isempty(f) || ~isempty(h)
        j = T.Tri(i, :);
        cord = T.Nodes(j, :);
        A = (cord(2:3,:) - repmat(cord(1,:),2,1))';
        Pcal = Pcord*A' + repmat(cord(1,:),size(Pcord,1),1);
        Area = abs(det(A));
        jp = Fd.NodePtrs( T.(U).TC(i, :) );
    end
    fval = f(Pcal(:, 1), Pcal(:, 2));
    I = Area*V'*(w'.*fval)/2;
    if ~isempty(f)
        for s = 1:ns
            if jp(s) > 0
                F(jp(s)) = F(jp(s)) + I(s);
                %F(jp(s)) = F(jp(s)) + w*(fval.*(P*C(:, s)))/2*Area;
            end
        end
    end
    
    if ~isempty(h)
        for s = 1:ns
            if jp(s) > 0
                F(jp(s)) = F(jp(s));
            end
        end
    end
end
end