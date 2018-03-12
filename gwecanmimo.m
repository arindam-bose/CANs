function [X, ISL, T] = gwecanmimo(N, M, alpha, beta, X0)
% X = gwecanmimo(N, M, alpha, beta) or X = gwecanmimo(N, M, alpha, beta, X0), G-WeCAN MIMO
%   N: length of the sequence X(i)
%   M: number of the sequence X(i)
%   alpha: N x 1, corresponding to weights w_k = gamma_k^2
%   beta: N x 1, corresponding to weights w_k = gamma_k^2
%   x0: N x 1, the initialization sequence
%   x: N x 1, the generated sequence

if nargin == 5
    X = X0;
else
    X = exp(1i * 2*pi * rand(N, M));
end
XPre = zeros(N, M);
iterDiff = norm(X - XPre);

alpha(1) = 0;
Alpha = toeplitz(alpha);
eigvalues = eig(Alpha);
alpha(1) = abs(min(eigvalues));
Alpha = toeplitz(alpha)/alpha(1);

beta(1) = 0;
Beta = toeplitz(beta);
eigvalues = eig(Beta);
beta(1) = abs(min(eigvalues));
Beta = toeplitz(beta)/beta(1);

kappa = alpha(1)*N*norm(Alpha,'fro')/(norm(Alpha,'fro')^2 + norm(Beta,'fro'));
k = 0;

figure(1);
figure(2);
tic;
while (iterDiff > 1e-3) && (k < 1000)
    k = k + 1;
    XPre = X;
    
    % step 2
    Xtilde = [X; zeros(N, M)]; % 2N x M
    D =  fft(Xtilde, 2*N, 1); % 2N x M
    Vp = zeros(2*N, M);
    for i = 1:2*N
        Vp(i,:) = sqrt(kappa) * D(i,:)/norm(D(i,:)); % 2N x M
    end
    
    % step 1
    G = ifft(Vp, 2*N, 1); % 2N x M
    X = exp(1i * angle(G(1:N,:))); % N x M
    % X = Gp(1:N,:); % if unimodularity is not a constraint
    % X = exp(1i * round(angle(Gp(1:N,:))/(2*pi/L))*(2*pi/L)); % if phase quantization to L levels, 
    % e.g. L = 256
    
    % iterDiff
    iterDiff = norm(X - XPre);
    
    % plots
    crr = xcorr2(X);
    figure(1);
    surf(-M+1:M-1, -N+1:N-1, 20*log(abs(crr)/N)/log(10));
    view([128, 17]);
    xlabel('index k');
    ylabel('index l');
    zlabel('crosscorrelation level (dB)');
    title({['G-WeCAN (MIMO): N = ' num2str(N) ', M = ' num2str(M)], ['Elapsed time: ' num2str(toc) ' sec']});
    
    ccrr = xcorr2(X,conj(X));
    figure(2);
    surf(-M+1:M-1, -N+1:N-1, 20*log(abs(ccrr)/N)/log(10));
    view([128, 17]);
    xlabel('index k');
    ylabel('index l');
    zlabel('complementary crosscorrelation level (dB)');
    title({['G-WeCAN (MIMO): N = ' num2str(N) ', M = ' num2str(M)], ['Elapsed time: ' num2str(toc) ' sec']});
    
    ISL = 20*log((norm(crr)^2-N^2)/sqrt(M*N^2))/log(10);
    T = toc;
    pause(0.000000000000000001);
end
end