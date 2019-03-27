function M = Bilinear1d(T, EdgeInd, U1, type1, U2, type2)
Quads1d = LoadQuad1d(4);
Tr1 = RefInfo1d(U1.Degree, type1, Quads1d);
Tr2 = RefInfo1d(U2.Degree, type2, Quads1d);

n1 = U1.Degree+1;
n2 = U2.Degree+1;
[edgei, edgej] = find(T.EgFlag(T.TrEg) == EdgeInd);
N = length(edgei);
indi = ones(n1*n2*N, 1);
indj = ones(n1*n2*N, 1);
valuek = zeros(n1*n2*N, 1);
indk = 1;
for i = 1:N
    iind = [edgej(i), edgej(i)+1-3*(edgej(i)==3), 3+edgej(i)];
    if U1.Degree == 1; ind1 = T.Tri(edgei(i), iind(1:2)); end
    if U1.Degree == 2; ind1 = [T.Tri(edgei(i), iind(1:2)), T.N+T.TrEg(edgei(i), edgej(i))]; end
    if U2.Degree == 1; ind2 = T.Tri(edgei(i), iind(1:2)); end
    if U2.Degree == 2; ind2 = [T.Tri(edgei(i), iind(1:2)), T.N+T.TrEg(edgei(i), edgej(i))]; end
    ind = T.Tri(edgei(i), :);
    cord = T.Node(ind, :);
    A = (cord(2:3, :) - repmat(cord(1, :),2,1))';   
    J = inv(A);
    cord = T.Node(ind(iind(1:2)), :);
    len = norm(cord(1, :)-cord(2, :), 2);
    for s = 1:n1
        for r = 1:n2
            indi(indk, 1) = ind1(s);
            indj(indk, 1) = ind2(r);
            switch type1
                case "mass"; I1 = Tr1.Ical(:, iind(s), edgej(i));
                case "dx"; I1 = ([Tr1.Ixcal(:, iind(s), edgej(i)), Tr1.Iycal(:, iind(s), edgej(i))]*J(:, 1));
                case "dy"; I1 = ([Tr1.Ixcal(:, iind(s), edgej(i)), Tr1.Iycal(:, iind(s), edgej(i))]*J(:, 2));
            end
            switch type2
                case "mass"; I2 = Tr2.Ical(:, iind(r), edgej(i));
                case "dx"; I2 = ([Tr2.Ixcal(:, iind(r), edgej(i)), Tr2.Iycal(:, iind(r), edgej(i))]*J(:, 1));
                case "dy"; I2 = ([Tr2.Ixcal(:, iind(r), edgej(i)), Tr2.Iycal(:, iind(r), edgej(i))]*J(:, 2));
            end
            valuek(indk, 1) = len/2*Quads1d.w*(I1.*I2);
            indk = indk+1;
        end
    end
end
M = sparse(indi, indj, valuek, U1.N, U2.N);
end