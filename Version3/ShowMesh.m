function ShowMesh(T, flag, color)
% This function show the mesh information
% flag controls the output
% flag(1) 0/1 contorl the info of nodes
% flag(2) 0/1 contorl the info of edges
% flag(3) 0/1 contorl the info of elements
if nargin < 3; color = 'kbrg'; end
if nargin < 2; flag = [0,0,0]; end

% Here plots the mesh
cordx = zeros(T.Ne, 2); cordy = zeros(T.Ne, 2);
cordx(:, 1) = T.Node(T.Edge(:, 1), 1);
cordx(:, 2) = T.Node(T.Edge(:, 2), 1);
cordy(:, 1) = T.Node(T.Edge(:, 1), 2);
cordy(:, 2) = T.Node(T.Edge(:, 2), 2);
plot(cordx', cordy', [color(1) '-']);

box off; set(gca, 'XTick', [], 'YTick', []);

% Show the info accroding to the flag
hold on
if flag(1)
   for i = 1:T.N
       text(T.Node(i,1), T.Node(i,2), num2str(i), 'color', color(2));
   end
end

if flag(2)
    for i = 1:T.Ne
        cord = mean(T.Node(T.Edge(i,:), :));
        text(cord(1), cord(2), num2str(i), 'color', color(3));
    end
end

if flag(3)
    for i = 1:T.Nt
        cord = mean(T.Node(T.Tri(i,:), :));
        text(cord(1), cord(2), num2str(i), 'color', color(4));
    end
end
hold off
axis(T.Shape);

end