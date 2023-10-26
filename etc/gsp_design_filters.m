function [hLP, hHP, hfLP, hfHP] = gsp_design_filters(lambda, L, showFlag, fig)
%gsp_design_filters
%   [hLP, hHP, hfLP, hfHP] = gsp_design_filters(lambda, L)
%   [hLP, hHP, hfLP, hfHP] = gsp_design_filters(lambda, L, showFlag)

if nargin < 3
  showFlag = 0;
end

if L > length(lambda)
  error('The order of the filters must be no larger than the number of eigenvalues.')
end

if showFlag && ~exist('fig', 'var')
  fig = figure;
end

% Frequency responses.
% Low-pass filter.
hfLP = real(sqrt(1-(regular(lambda*(2/max(lambda))).^2)));
% High-pass filter.
hfHP = regular(lambda*(2/max(lambda)));

Psi = fliplr(vander(lambda));
PsiInv = pinv(Psi(:, 1:L));

hLP = PsiInv*hfLP;
hHP = PsiInv*hfHP;

if showFlag
  figure(fig)
  clf

  subplot(211)
  hold on
  plot(lambda, hfLP, 'o-')
  plot(lambda, hfHP, 'o-')
  plot(lambda, Psi(:, 1:L)*hLP, 'x-')
  plot(lambda, Psi(:, 1:L)*hHP, 'x-')
  hold off
  legend('True LP', 'True HP', 'Approx. LP', 'Approx. HP')

  subplot(212)
  hold on
  plot(lambda, hfLP.*hfHP, 'o-')
  plot(lambda, (Psi(:, 1:L)*hLP).*(Psi(:, 1:L)*hHP), 'x-')
  hold off
  legend('True LP-HP Hadamard', 'Approx. LP-HP Hadamard')

  drawnow
end

end

function y = regular(val)
% Copied from gsp_design_regular.m in the GSP toolbox https://lts2.epfl.ch/gsp/.

% Filter degree.
d = 4;

if d == 0
  y = sin(pi/4*val);
else
  output = sin(pi*(val-1)/2);
  for k = 2:d
    output = sin(pi*output/2);
  end
  y = sin(pi/4*(1+output));
end

end
