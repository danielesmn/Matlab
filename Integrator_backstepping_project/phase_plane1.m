function [y] = phase_plane1(x,color)

x1 = x(1, :);
x2 = x(2, :);
x3 = x(3, :);

plot3(x1, x2, x3, color, 'LineWidth', 1.5);

grid on
hold on

u = gradient(x1);
v = gradient(x2);
w = gradient(x3);
quiver3(x1, x2, x3, u, v, w);

end

