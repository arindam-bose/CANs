function [x, ISL] = wecansiso(N, gamma, x0)
% x = wecansiso(N, gamma) or x = wecansiso(N, gamma, x0), WeCAN SISO
%   N: length of the sequence x
%   gamma: N x 1, corresponding to weights w_k = gamma_k^2
%   x0: N x 1, the initialization sequence
%   x: N x 1, the generated sequence

if nargin == 3
    x = x0;
else
    x = exp(1i * 2*pi * rand(N, 1));
end
xPre = zeros(N,1);
iterDiff = norm(x - xPre);

gamma(1) = 0;
Gamma = toeplitz(gamma);
eigvalues = eig(Gamma);
gamma(1) = abs(min(eigvalues));
Gamma = toeplitz(gamma)/gamma(1);
C = sqrtm(Gamma);

Alpha = zeros(2*N, N); % Alpha(p,:) is alpha_p

figure;
tic;
while(iterDiff > 1e-4)
    xPre = x;
    
    % step 1
    Z = [C.' .* kron(x,ones(1,N)); zeros(N,N)]; % 2N x N
    F = fft(Z); % 2N*N, p(th) row corresponds to alpha_p
    for p = 1:(2*N)
        Alpha(p,:) = sqrt(N) * F(p,:) / norm(F(p,:));
    end
    
    % step 2
    Nu = ifft(Alpha); % 2N x N
    MuNu = C' .* Nu(1:N,:); % N x N
    for n = 1:N
        x(n) = exp(1i * angle(sum(MuNu(n,:))));
    end
    
    % stop criterion
    iterDiff = norm(x - xPre);
    
    % plot   
    crr = xcorr(x);
    % surf(-M+1:M-1, -N+1:N-1, 20*log(abs(crr)));
    plot(-N+1:N-1, 20*log(abs(crr)/N)/log(10),'LineWidth',1.3);
    xlabel('index k');
    ylabel('autocorrelation level (dB)');
    xlim([-N+1 N-1]);
    title({['WeCAN (SISO): N = ' num2str(N) ', M = 1'], ['Elapsed time: ' num2str(toc) ' sec']});
    ISL = norm(crr)^2-N^2;
    pause(0.000000000000000001);
end