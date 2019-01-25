function ShowMesh(T, flag)
    if nargin < 2; flag = [0 0 0];end
    triplot(T.Tri, T.Nodes(:,1), T.Nodes(:,2));
    box off; set(gca, 'XTick', [], 'YTick', []);
    
    if flag(1)
        for i = 1:T.Nt
            cord = T.Nodes(T.Tri(i,:),:);
            text(sum(cord(:,1))/3, sum(cord(:,2))/3, num2str(i), 'Color', 'k');
        end
    end
    if flag(2)
        for i = 1:T.N
            text(T.Nodes(i,1), T.Nodes(i,2), num2str(i), 'Color', 'r');
        end
    end
    if flag(3)
        for i = 1:T.Ne
            cord = T.Nodes(T.Edge(i,:),:);
            text(sum(cord(:,1))/2, sum(cord(:,2))/2, num2str(i), 'Color', 'g');
        end
    end
end