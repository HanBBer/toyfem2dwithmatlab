function K = FEMatrix(T, Fd1, d1, Fd2, d2)
U1 = Fd1.Space;
if nargin < 4; U2 = U1; Fd2 = Fd1; d2 = d1; 
else; U2 = Fd2.Space;
end
[w, P, Px, Py, C01, C02] = LoadQuad();
if (T.(U1).Property == "P1");   C1 = C01;
else;  C1 = C02; end
if (T.(U1).Property == "P1");   C2 = C01;
else;  C2 = C02; end
Nf1 = size(Fd1.FNodePtrs, 1);
Nf2 = size(Fd2.FNodePtrs, 1);
K = zeros(Nf1, Nf2);
ns = size(T.(U1).TC, 2); nr = size(T.(U2).TC, 2);
for i = 1:T.Nt
    j = T.Tri(i, :);
    cord = T.Nodes(j,:);
    A = (cord(2:3,:) - repmat(cord(1,:),2,1))';
    J = inv(A);
    Area = abs(det(A));
    j1 = T.(U1).TC(i, :); j2 = T.(U2).TC(i, :);
    jp1 = Fd1.NodePtrs(j1); jp2 = Fd2.NodePtrs(j2);
    for s = 1:ns
        if jp1(s) > 0
            c1 = C1(:, s);
            for r = 1:nr
                if jp2(r) > 0
                    c2 = C2(:, r);
                    if d1 == "mass"
                        I1 = P*c1;
                    else
                        V1 = [Px*c1 Py*c1];
                        switch d1
                            case {"dx"}; I1 = V1*J(:,1);
                            case {"dy"}; I1 = V1*J(:,2);
                            case {"nabla"}; I1 = V1/A;
                        end
                    end
                    if d2 == "mass"
                        I2 = P*c2;
                    else
                        V2 = [Px*c2 Py*c2];
                        switch d2
                            case {"dx"}; I2 = V2*J(:,1);
                            case {"dy"}; I2 = V2*J(:,2);
                            case {"nabla"}; I2 = V2/A;
                        end
                    end

                    I = sum(w*(I1.*I2)/2);
                    K(jp1(s), jp2(r)) = K(jp1(s), jp2(r)) + I*Area; 
                end
            end
        end
    end
end

end