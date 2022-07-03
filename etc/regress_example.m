clearvars

xmin = 1.4e5;
xmax = 2e5;

xp = xmin + (xmax-xmin)*rand(2,1);
xp = sort(xp);

y0 = 1.1;
y1 = 0.3;

m = (y1-y0) / (xp(2)-xp(1));
b = 8;

N = 100;

x = xmin + (xmax-xmin)*rand(N,1);
y = m*x + b;

x = x + xmin/10*randn(N, 1);

scatter(x, y)

coeffs = regress(y, [ones(N, 1) x]);

x_axis = linspace(xmin, xmax);
sol = coeffs(1) + coeffs(2)*x_axis;

hold on
plot(x_axis, sol)
hold off

grid on
box on
title(sprintf('%.2d %.2d', coeffs))
xlabel('x')
ylabel('y')
