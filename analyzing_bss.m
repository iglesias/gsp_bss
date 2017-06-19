clearvars

% Number of graph nodes.
N = 3;
% Number of filter coefficients.
L = 1;

x1 = sym('x1', [N 1], 'real');
x2 = sym('x2', [N 1], 'real');

h1 = sym('h1', [L 1]', 'real');
h2 = sym('h2', [L 1]', 'real');

Z1 = x1*h1';
Z2 = x2*h2';

Z = Z1+Z2;

vec(Z)
