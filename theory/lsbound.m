% Modified 31/01/2019

N = 50;
L = N;

r = 2;
K = 3;

alpha = 1:0.1:3;
gamma = sqrt(N*log(N*L/2) + alpha*log(L));

mumax2 = L/K;
muh2 = 16/9*mumax2*K;

bound_wo_max = r^2*(log(L))^2*log(gamma);
bound = bound_wo_max*max(mumax2*K, muh2*N);

subplot(211)
plot(alpha, 1-L.^(-alpha+1))
xlabel('alpha'), ylabel('Probability 1-L.\^(-alpha+1)')
grid on
subplot(212)
plot(1-L.^(-alpha+1), bound_wo_max)
ylabel('bound'), xlabel('Probability 1-L.\^(-alpha+1)')
grid on


% See Blind Deconvolution Meets Blind Demixing: Algorithms and Performance Bounds
