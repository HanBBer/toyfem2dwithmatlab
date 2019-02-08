function Quads = LoadQuad1d(degree)
% This function provides the 1-d quadrature info on ]0,1[
% including the weights w and corres

if (degree == 1)||(degree == 0)
    Quads.w = 2;
    Quads.px = 0;
end

if degree == 2
    Quads.w = [1, 1];
    Quads.px = [0.577350269189626; -0.577350269189626];
end

if degree == 3
    Quads.w = [0.888888888888889, 0.555555555555556*ones(1, 2)];
    Quads.px = [0; 0.774596669241483*[1;-1]];
end

if degree == 4
    Quads.w = [0.652145154862546*ones(1, 2), 0.347854845137454*ones(1, 2)];
    Quads.px = [0.339981043584856*[1;-1]; 0.861136311594053*[1;-1]];
end

if degree == 5
    Quads.w = [0.568888888888889, 0.478628670499366*ones(1, 2), 0.236926885056189*ones(1, 2)];
    Quads.px = [0;
        0.059715871789770; 0.470142064105115; 0.470142064105115;
        0.797426985353087; 0.101286507323456; 0.101286507323456];
end

if degree == 6
    Quads.w = [0.116786275726379*ones(1, 3), 0.050844906370207*ones(1,3), 0.082851075618374*ones(1,6)];
    Quads.px = [0; 0.538469310105683*[1;-1]; 0.906179845938664*[1;-1]];
end


end