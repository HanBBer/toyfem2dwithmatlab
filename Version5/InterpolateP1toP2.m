function v = InterpolateP1toP2(T, u)
% This function transfer the P1 FE function to P2 FE function
v = [u; 0.5*( u(T.Edge(:,1)) + u(T.Edge(:,2)) )];
end