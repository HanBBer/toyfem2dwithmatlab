load('data.mat')
A = 3;


figure(1)
for iter = 6
    for i = 1:3
        for j = 1:3
            subplot(3,3, 3*(i-1)+j)
            result = RESULT{i, j};
            T0 = result.TF;
            T = result.data{iter,1};
            sol = result.data{iter,2};
            trisurf(T.Tri, T0.Node(:,1)+A*(T.Node(:, 1)-T0.Node(:,1)), T0.Node(:,2)+A*(T.Node(:, 2)-T0.Node(:,2)), sol(2*result.NFu+1:2*result.NFu+result.NFp),...
                'FaceColor', 'interp', 'EdgeColor', 'interp');
            axis([0, 6, -1, 1])
            axis equal
            box off; set(gca, 'XTick', [], 'YTick', [], 'xcolor', 'white', 'ycolor', 'white');
            view(2);
        end
    end
    %pause(0.05)
end


style = {'k-','b-.','r.'};
lengthf = 1000;
widthf = 400;
lengthp = 250;
intervalp = 70;
fgx = (lengthf-3*lengthp-2*intervalp)/lengthf/2;
fgy = (widthf-lengthp)/widthf/2;
pgx = lengthp/lengthf;
itvx = intervalp/lengthf;
pgy = lengthp/widthf;
%%
figure(2)
set (gcf,'position',[100 100 lengthf widthf]);
for j = 1:3
    subplot('position', [fgx+(pgx+itvx)*(j-1), fgy, pgx, pgy])
    hold on
    for i = 1:3
        result = RESULT{i, j};
        TS = result.TS;
        sol = result.data{4,2};
        d = sol(2*result.NFu+result.NFp+3*result.NSu+1:2*result.NFu+result.NFp+4*result.NSu);
        posSP1 = TS.Edge(TS.EgFlag == 1, :);
        posSP1 = unique(posSP1(:));
        [~, n] = sort(TS.Node(posSP1, 1));
        posSP1 = posSP1(n);
        plot(TS.Node(posSP1,1), d(posSP1), style{i})
        axis([0, 6, -0.1, 0.1])
    end
    box on
    legend({'h = 0.2', 'h = 0.1', 'h = 0.05'}, 'location', 'south');
    title(['\Deltat = ', num2str(0.2/2^j)]);
    hold off
end

%%
figure(3)
set (gcf,'position',[100 100 lengthf widthf]);
for i = 1:3
    subplot('position', [fgx+(pgx+itvx)*(i-1), fgy, pgx, pgy])
    hold on
    for j = 1:3
        result = RESULT{i, j};
        TS = result.TS;
        sol = result.data{4,2};
        d = sol(2*result.NFu+result.NFp+3*result.NSu+1:2*result.NFu+result.NFp+4*result.NSu);
        posSP1 = TS.Edge(TS.EgFlag == 1, :);
        posSP1 = unique(posSP1(:));
        [~, n] = sort(TS.Node(posSP1, 1));
        posSP1 = posSP1(n);
        plot(TS.Node(posSP1,1), d(posSP1), style{j})
        axis([0, 6, -0.1, 0.1])
    end
    box on
    legend({'\Deltat = 0.1', '\Deltat = 0.05', '\Deltat = 0.025'}, 'location', 'south');
    title(['h = ', num2str(0.4/2^i)]);
    hold off
end


%%
figure(4)
lengthf = 800;
widthf = 900;
lengthp = 600;
widthp = 100;
intervalp = 50;
fgx = (lengthf-lengthp)/lengthf/2;
fgy = (widthf-5*widthp-4*intervalp)/widthf/2;
pgx = lengthp/lengthf;
itvy = intervalp/widthf;
pgy = widthp/widthf;

result = RESULT{3, 3};
T0 = result.TF;
set (gcf,'position',[100 00 lengthf/1.5 widthf/1.5]);
for i = 1:5
axes('position', [fgx, fgy+(pgy+itvy)*(i-1), pgx, pgy])
T = result.data{7-i, 1};
sol = result.data{7-i, 2};
trisurf(T.Tri, T0.Node(:,1)+A*(T.Node(:, 1)-T0.Node(:,1)), T0.Node(:,2)+A*(T.Node(:, 2)-T0.Node(:,2)), sol(2*result.NFu+1:2*result.NFu+result.NFp),...
    'FaceColor', 'interp', 'EdgeColor', 'interp');
hold on
trisurf(T.Tri, T0.Node(:,1)-A*(T.Node(:, 1)-T0.Node(:,1)), -T0.Node(:,2)-A*(T.Node(:, 2)-T0.Node(:,2)), sol(2*result.NFu+1:2*result.NFu+result.NFp),...
    'FaceColor', 'interp', 'EdgeColor', 'interp');
hold off
axis([0, 6, -0.5, 0.5])
%box off; set(gca, 'XTick', [], 'YTick', [], 'xcolor', 'white', 'ycolor', 'white');
axis off
view(2);
title(['t = ', num2str(0.1*(6-i)), ' s'])
end


%%
figure(5)
lengthf = 800;
widthf = 900;
lengthp = 600;
widthp = 100;
intervalp = 50;
fgx = (lengthf-lengthp)/lengthf/2;
fgy = (widthf-5*widthp-4*intervalp)/widthf/2;
pgx = lengthp/lengthf;
itvy = intervalp/widthf;
pgy = widthp/widthf;

result = RESULT{1, 3};
T0 = result.TF;
set (gcf,'position',[100 00 lengthf/1.5 widthf/1.5]);
for i = 1:5
axes('position', [fgx, fgy+(pgy+itvy)*(i-1), pgx, pgy])
T = result.data{13-2*i, 1};
sol = result.data{13-2*i, 2};
triplot(T.Tri, T0.Node(:,1)+A*(T.Node(:, 1)-T0.Node(:,1)), T0.Node(:,2)+A*(T.Node(:, 2)-T0.Node(:,2)));
hold on
triplot(T.Tri, T0.Node(:,1)-A*(T.Node(:, 1)-T0.Node(:,1)), -T0.Node(:,2)-A*(T.Node(:, 2)-T0.Node(:,2)));
hold off
axis([0, 6, -0.5, 0.5])
%box off; set(gca, 'XTick', [], 'YTick', [], 'xcolor', 'white', 'ycolor', 'white');
axis off
view(2);
title(['t = ', num2str(0.1*(12-2*i)), ' s'])
end


load('data0.mat')
iter = 8;
TF0R = Reference.TF;
TS0R = Reference.TS;
solR = Reference.data{iter, 2};
uR = [solR(1:Reference.NFu), solR(Reference.NFu+1:2*Reference.NFu)];
etaR = [solR(2*Reference.NFu+Reference.NFp+2*Reference.NSu+(1:Reference.NSu)), solR(2*Reference.NFu+Reference.NFp+3*Reference.NSu+(1:Reference.NSu))];
for i = 1:3
    for j = 1:3
        result = RESULT{i, j};
        sol = result.data{iter, 2};
        u = [sol(1:result.NFp), sol(result.NFu+1:result.NFu+result.NFp)];
        eta = [sol(2*result.NFu+result.NFp+2*result.NSu+(1:result.NSp)), sol(2*result.NFu+result.NFp+3*result.NSu+(1:result.NSp))];
        
        TF0 = result.TF;
        ix = ismember(TF0R.Node(:,1), TF0.Node(:,1));
        iy = ismember(TF0R.Node(:,2), TF0.Node(:,2));
        ind = ix & iy;
        eF = uR(ind,:) - u;
        
        TS0 = result.TS;
        ix = ismember(TS0R.Node(:,1), TS0.Node(:,1));
        iy = ismember(TS0R.Node(:,2), TS0.Node(:,2));
        ind = ix & iy;
        eS = etaR(ind,:) - eta;
        
        TF = result.data{iter, 1};
        PF = P1Fespace(TF);
        EF = symBilinear(PF, 'mass', []) + symBilinear(PF, 'nabla', []);
        
        PS = P1Fespace(TS0);
        ES = symBilinear(PS, 'mass', []) + symBilinear(PS, 'nabla', []);
        
        fprintf('Error of u: \t%.4f\n', norm(EF*eF(:,1))+norm(EF*eF(:,2)))
        fprintf('Error of eta: \t%4e\n', norm(ES*eS(:,1))+norm(ES*eS(:,2)))
    end
end