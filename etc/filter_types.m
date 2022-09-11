clearvars

L = 3;
R = 2;

figure
hold on

for i = 1:100
  % case DataDistribution.Uniform
%   h = rand(L, R);
  % case DataDistribution.Normal
%   h = randn(L, R);
% case DataDistribution.HeatKernel
  h = exp(-(0.3*rand(R, 1)+0.7)*(1:L))';
  % Normalize filter taps.
  h = h ./ repmat(norms(h, 2, 1), L, 1);
  plot(h, 'o-')
end
