function [X, ISL] = gcanmimo(N, M, X0)
% X = gcanmimo(N, M) or X = gcanmimo(N, M, X0), G-CAN MIMO
%   N: length of the sequence X(i)
%   M: number of the sequence X(i)
%   X0: N x M, the initialization set of sequences
%   X: N x M, the generated set of sequences

if nargin == 3
    X = X0;
else
    X = exp(1i * 2*pi * rand(N, M)); % N x M
end

XPre = zeros(N, M);
iterDiff = norm(X - XPre);
k = 0;

figure;
tic;
while (iterDiff>1e-6) && (k<1000)
    k = k + 1;
    XPre = X;
    
    % step 2
    Xtilde = [X; zeros(N, M)]; % 2N x M
    D =  1/sqrt(N) * fft(Xtilde, 2*N, 1); % 2N x M
    V = zeros(2*N, M);
    for i = 1:2*N
        V(i,:) = sqrt(1/2) * D(i,:)/norm(D(i,:)); % 2N x M
    end
    
    % step 1
    G = sqrt(N) * ifft(V, 2*N, 1); % 2N x M
    X = exp(1i * angle(G(1:N,:))); % N x M
    % X = Gp(1:N,:); % if unimodularity is not a constraint
    % X = exp(1i * round(angle(Gp(1:N,:))/(2*pi/L))*(2*pi/L)); % if phase quantization to L levels, 
    % e.g. L = 256
    
    % iterDiff
    iterDiff = norm(X - XPre);
    
    % plots
    crr = xcorr2(X);
    surf(-M+1:M-1, -N+1:N-1, 20*log(abs(crr)/N)/log(10));
    xlabel('index k');
    ylabel('index l');
    zlabel('autocorrelation level (dB)');
    title({['G-CAN (MIMO): N = ' num2str(N) ', M = ' num2str(M)], ['Elapsed time: ' num2str(toc) ' sec']});
    ISL = norm(crr)^2-N^2;
    pause(0.000000000000000001);
end
end