% this script demonstrates the function 'mutual_coherence'
% Written by Dr. Yoash Levron, Technion, Israel, 2015

clc;

% first example:
M=20;  N = 500;
A = 10*randn(M,N) + 10i*randn(M,N);
u = mutual_coherence(A)
lower_bound_on_mutual_coherence = ((N-M)/(M*(N-1)))^0.5

% second example:
M=2;  N = 4;
A = randn(M,N);
u = mutual_coherence(A)
lower_bound_on_mutual_coherence_for_2D_columns = cos(pi/N)

% Plot the normalized columns of A (for the second example):
B = bsxfun(@rdivide,A,sqrt(sum(A.*conj(A),1))); % normalize the columns
plot(B(1,:),B(2,:),'kx','linewidth',7);
xlim([-1 1]); ylim([-1 1]);
title('normalized columns of A');
