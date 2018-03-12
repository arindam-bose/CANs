function [X, ISL, T] = wecanmimo(N, M, gamma, X0)
% X = wecanmimo(N, M, gamma) or X = wecanmimo(N, M, gamma, X0), WeCAN MIMO
%   N: length of the sequence X(i)
%   M: number of the sequence X(i)
%   gamma: N x 1, corresponding to weights w_k = gamma_k^2
%   x0: N x 1, the initialization sequence
%   x: N x 1, the generated sequence

if nargin == 4
    X = X0;
else
    X = exp(1i * 2*pi * rand(N, M));
end
XPre = zeros(N, M);
iterDiff = norm(X - XPre);

gamma(1) = 0;
Gamma = toeplitz(gamma);
eigvalues = eig(Gamma);
gamma(1) = abs(min(eigvalues));
Gamma = toeplitz(gamma)/gamma(1);
C = sqrtm(Gamma);
U = zeros(2*N, M*N); % Alpha(p,:) is alpha_p
k = 0;

figure;
tic;
while(iterDiff > 1e-3)  %&& (k < 5000)
    k = k + 1;
    XPre = X;
    
    % step 1
    for i = 1:M
        XiTilde = [(C.').*(diag(X(:,i))*ones(N,N)); zeros(N,N)];
        if i == 1
            FTilde = XiTilde;
        else
            FTilde = [FTilde XiTilde]; % 2N x MN
        end
    end
    F = sqrt(2*N) * fft(FTilde); % 2N x MN
    for p = 1:(2*N)
        fp = reshape(F(p,:), N, M);
        [U1,~,U2] = svd(fp', 'econ');
        U(p,:) = reshape(U2*U1', M*N, 1); % 2N x MN
    end
    
    % step 2
    G = sqrt(gamma(1)*N) .* U; % 2N x MN
    Nu = 1/sqrt(2*N) * ifft(G); % 2N x MN
    for m = 1:M
        MuNu = C' .* Nu(1:N,(m-1)*N+1:m*N); % N x N
        for n = 1:N
            X(n,m) = exp(1i * angle(sum(MuNu(n,:))));
        end
    end
    
    % stop criterion
    iterDiff = norm(X - XPre);
    
    % plot   
    crr = xcorr2(X);
    surf(-M+1:M-1, -N+1:N-1, 20*log(abs(crr)/N)/log(10));
    view([128, 17]);
    xlabel('index k');
    ylabel('index l');
    zlabel('autocorrelation level (dB)');
    title({['WeCAN (MIMO): N = ' num2str(N) ', M = ' num2str(M)], ['Elapsed time: ' num2str(toc) ' sec']});
    ISL = 20*log((norm(crr)^2-N^2)/sqrt(M*N^2))/log(10);
    T = toc;
    pause(0.000000000000000001);
end