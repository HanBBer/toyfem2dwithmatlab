function Tr = RefInfo(degree, Quads2d)
% This function prepare the Quad info on reffence Trigular
% accroding to the degree of FEspace
if degree == 1
    Tr.P = [ones(size(Quads2d.px)), Quads2d.px, Quads2d.py];
    Tr.Px = [zeros(size(Quads2d.px)), ones(size(Quads2d.px)), zeros(size(Quads2d.px))];
    Tr.Py = [zeros(size(Quads2d.px)), zeros(size(Quads2d.px)), ones(size(Quads2d.px))];
    Rcord = [0, 0; 1, 0; 0, 1];
    Tr.C = [ones(size(Rcord,1),1), Rcord];   
end
if degree == 2
    Tr.P = [ones(size(Quads2d.px)), Quads2d.px, Quads2d.py,...
        Quads2d.px.^2, Quads2d.px.*Quads2d.py, Quads2d.py.^2];
    Tr.Px = [zeros(size(Quads2d.px)), ones(size(Quads2d.px)), zeros(size(Quads2d.px)),...
        2*Quads2d.px, Quads2d.py, zeros(size(Quads2d.px))];
    Tr.Py = [zeros(size(Quads2d.px)), zeros(size(Quads2d.px)), ones(size(Quads2d.px)),...
        zeros(size(Quads2d.px)), Quads2d.px, 2*Quads2d.py];
    Rcord = [0, 0; 1, 0; 0, 1; 1/2, 0; 1/2, 1/2; 0, 1/2];
    Tr.C = [ones(size(Rcord,1),1), Rcord, Rcord(:,1).^2, Rcord(:,1).*Rcord(:,2), Rcord(:,2).^2];
end

end