function [X, Y, Z] = cube2line(cube)
vetex = cube.d3Vertices(1, :);
X = [vetex(1),  vetex(2), vetex(3), vetex(4), vetex(1), vetex(5), vetex(6), vetex(7), vetex(8), vetex(5), ...
    vetex(6), vetex(2), vetex(3), vetex(7), vetex(8), vetex(4)];
vetex = cube.d3Vertices(2, :);
Y = [vetex(1),  vetex(2), vetex(3), vetex(4), vetex(1), vetex(5), vetex(6), vetex(7), vetex(8), vetex(5), ...
    vetex(6), vetex(2), vetex(3), vetex(7), vetex(8), vetex(4)];
vetex = cube.d3Vertices(3, :);
Z = [vetex(1),  vetex(2), vetex(3), vetex(4), vetex(1), vetex(5), vetex(6), vetex(7), vetex(8), vetex(5), ...
    vetex(6), vetex(2), vetex(3), vetex(7), vetex(8), vetex(4)];
end