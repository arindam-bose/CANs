function [x, ISL] = cansiso(N, x0)
% x = cansiso(N) or x = cansiso(N, x0), CAN SISO
%   N: length of the sequence x
%   x0: N x 1, the initialization sequence
%   x: N x 1, the generated sequence

if nargin == 2
    x = x0;
else
    x = exp(1i * 2*pi * rand(N,1));
end

xPre = zeros(N, 1);
iterDiff = norm(x - xPre);
k = 0;

figure;
tic;
while (iterDiff>1e-6) && (k<1000)
    k = k + 1;
    xPre = x;
    
    % step 2
    xTilde = [x; zeros(N, 1)]; % 2N x 1
    d = 1/sqrt(2*N) * fft(xTilde); % 2N x 1
    v = sqrt(1/2) * exp(1i * angle(d)); % 2N x 1
    
    % step 1
    g = sqrt(2*N) * ifft(v); % 2N x 1    
    x = exp(1i * angle(g(1:N))); % N x 1
    % x = g(1:N); % if unimodularity is not a constraint
    % x = exp(1i * round(angle(g(1:N))/(2*pi/L))*(2*pi/L)); % if phase quantization to L levels,
    % e.g. L = 256
    
    % iterDiff
    iterDiff = norm(x - xPre);

    % plot 
    crr = xcorr(x);
    % surf(-M+1:M-1, -N+1:N-1, 20*log(abs(crr)));
    plot(-N+1:N-1, 20*log(abs(crr)/N)/log(10),'LineWidth',1.3);
    xlabel('index k');
    ylabel('autocorrelation level (dB)');
    xlim([-N+1 N-1]);
    title({['CAN (SISO): N = ' num2str(N) ', M = 1'], ['Elapsed time: ' num2str(toc) ' sec']});
    ISL = norm(crr)^2-N^2;
    pause(0.000000000000000001);
end