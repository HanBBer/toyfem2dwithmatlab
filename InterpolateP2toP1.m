function v = InterpolateP2toP1(T, u)
% This function transfer the P2 FE function to P1 FE function
v = u(1:T.N);
end