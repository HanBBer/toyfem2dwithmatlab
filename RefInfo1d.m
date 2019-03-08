function Tr = RefInfo1d(degree, type, Quads1d)
% This function returns the info for 1-d integral on Mesh
nquads1d = size(Quads1d.px, 1);
if type == "mass"
    Tr.Ical = zeros(nquads1d, 3*degree, 3);
else
    Tr.Ixcal = zeros(nquads1d, 3*degree, 3);
    Tr.Iycal = zeros(nquads1d, 3*degree, 3);
end
if degree == 1
    Rcord = [0, 0; 1, 0; 0, 1];
    C = [ones(size(Rcord,1),1), Rcord]; 
    if type == "mass"
        Tr.Ical(:, :, 1) = [ones(nquads1d, 1), 1/2+Quads1d.px/2, zeros(nquads1d, 1)]/C;
        Tr.Ical(:, :, 2) = [ones(nquads1d, 1), 1/2-Quads1d.px/2, 1/2+Quads1d.px/2]/C;
        Tr.Ical(:, :, 3) = [ones(nquads1d, 1), zeros(nquads1d, 1), 1/2-Quads1d.px/2]/C;
    else
        Tr.Ixcal(:, :, 1) = [zeros(nquads1d, 1), ones(nquads1d, 1), zeros(nquads1d, 1)]/C;
        Tr.Ixcal(:, :, 2) = [zeros(nquads1d, 1), ones(nquads1d, 1), zeros(nquads1d, 1)]/C;
        Tr.Ixcal(:, :, 3) = [zeros(nquads1d, 1), ones(nquads1d, 1), zeros(nquads1d, 1)]/C;
        Tr.Iycal(:, :, 1) = [zeros(nquads1d, 1), zeros(nquads1d, 1), ones(nquads1d, 1)]/C;
        Tr.Iycal(:, :, 2) = [zeros(nquads1d, 1), zeros(nquads1d, 1), ones(nquads1d, 1)]/C;
        Tr.Iycal(:, :, 3) = [zeros(nquads1d, 1), zeros(nquads1d, 1), ones(nquads1d, 1)]/C;
    end
end
if degree == 2
    Rcord = [0, 0; 1, 0; 0, 1; 1/2, 0; 1/2, 1/2; 0, 1/2];
    C = [ones(size(Rcord,1),1), Rcord, Rcord(:,1).^2, Rcord(:,1).*Rcord(:,2), Rcord(:,2).^2];
    if type == "mass"
        Tr.Ical(:, :, 1) = [ones(nquads1d, 1), 1/2+Quads1d.px/2, zeros(nquads1d, 1), (1/2+Quads1d.px/2).^2, zeros(nquads1d, 2)]/C;
        Tr.Ical(:, :, 2) = [ones(nquads1d, 1), 1/2-Quads1d.px/2, 1/2+Quads1d.px/2, (1/2-Quads1d.px/2).^2, (1/2-Quads1d.px/2).*(1/2+Quads1d.px/2), (1/2+Quads1d.px/2).^2]/C;
        Tr.Ical(:, :, 3) = [ones(nquads1d, 1), zeros(nquads1d, 1), 1/2-Quads1d.px/2, zeros(nquads1d, 2), (1/2-Quads1d.px/2).^2]/C;
    else
        Tr.Ixcal(:, :, 1) = [zeros(nquads1d, 1), ones(nquads1d, 1), zeros(nquads1d, 1), 1+Quads1d.px, zeros(nquads1d, 2)]/C;
        Tr.Ixcal(:, :, 2) = [zeros(nquads1d, 1), ones(nquads1d, 1), zeros(nquads1d, 1), 1-Quads1d.px, 1/2+Quads1d.px/2, zeros(nquads1d, 1)]/C;
        Tr.Ixcal(:, :, 3) = [zeros(nquads1d, 1), ones(nquads1d, 1), zeros(nquads1d, 2), 1/2-Quads1d.px/2, zeros(nquads1d, 1)]/C;
        Tr.Iycal(:, :, 1) = [zeros(nquads1d, 1), zeros(nquads1d, 1), ones(nquads1d, 1), zeros(nquads1d, 1), 1/2+Quads1d.px/2, zeros(nquads1d, 1)]/C;
        Tr.Iycal(:, :, 2) = [zeros(nquads1d, 1), zeros(nquads1d, 1), ones(nquads1d, 1), zeros(nquads1d, 1), 1/2-Quads1d.px/2, 1+Quads1d.px]/C;
        Tr.Iycal(:, :, 3) = [zeros(nquads1d, 1), zeros(nquads1d, 1), ones(nquads1d, 1), zeros(nquads1d, 2), 1-Quads1d.px]/C;
    end
end

end