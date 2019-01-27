function [w, P, Px, Py, C1, C2]=LoadQuad()
% This function provides the data used in calculation
w = [0.116786275726379*ones(1, 3), 0.050844906370207*ones(1,3), 0.082851075618374*ones(1,6)];
px = [0.501426509658179; 0.249286745170910; 0.249286745170910;
    0.873821971016996; 0.063089014491502; 0.063089014491502;
    0.053145049844817; 0.053145049844817; 0.310352451033784;
    0.310352451033784; 0.636502499121399; 0.636502499121399];
py = [0.249286745170910; 0.249286745170910; 0.501426509658179;
    0.063089014491502; 0.063089014491502; 0.873821971016996;
    0.310352451033784; 0.636502499121399; 0.636502499121399;
    0.053145049844817; 0.053145049844817; 0.310352451033784;];
P = [ones(12,1), px, py, px.^2, px.*py, py.^2];
Px = [zeros(12,1), ones(12,1), zeros(12,1), 2*px, py, zeros(12,1)];
Py = [zeros(12,1), zeros(12,1), ones(12,1), zeros(12,1), px, 2*py];

cord = [0, 0; 1, 0; 0, 1];
M = [ones(3,1), cord];
C1 = zeros(6, 6);
C1(1:3, 1:3) = inv(M); % Coefficients of basis funcitons
cord = [0, 0; 1, 0; 0, 1; 1/2, 0; 1/2, 1/2; 0, 1/2]; % The nodal points
M = [ones(6,1), cord, cord(:, 1).^2, cord(:, 1).*cord(:, 2), cord(:, 2).^2];
C2 = inv(M); % Coefficients of basis funcitons


end