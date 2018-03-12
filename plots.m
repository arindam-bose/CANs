close all;
clear all;
clc;

%%
M = 3;
Ns = [10, 30, 100, 300, 1000];
P = 0.3;
ISL = zeros(length(Ns), 3);
time = zeros(length(Ns), 3);

for i = 1:length(Ns)
    N = Ns(i);
    alpha = [ones(N*P,1); zeros(N*(1-P),1)];
    beta = [ones(N*P,1); zeros(N*(1-P),1)];
    X0 = exp(1i * 2*pi * rand(N, M));
    [~, ISL(i, 1), time(i, 1)] = canmimo(N, M, X0);
    disp(['CAN    : N = ' num2str(N) ', M = ' num2str(M) ', ISL = ' num2str(ISL(i, 1)) ' dB, Time = ' num2str(time(i, 1)) ' sec']);
    [~, ISL(i, 2), time(i, 2)] = wecanmimo(N, M, alpha, X0);
    disp(['WeCAN  : N = ' num2str(N) ', M = ' num2str(M) ', ISL = ' num2str(ISL(i, 2)) ' dB, Time = ' num2str(time(i, 2)) ' sec']);
    [~, ISL(i, 3), time(i, 3)] = gwecanmimo(N, M, alpha, beta, X0);
    disp(['G-WeCAN: N = ' num2str(N) ', M = ' num2str(M) ', ISL = ' num2str(ISL(i, 3)) ' dB, Time = ' num2str(time(i, 3)) ' sec']);
    disp('-----------------------------------------------------------------');
end

figure;
semilogx(Ns, ISL(:,1), '--o', 'LineWidth', 2);
hold on;
semilogx(Ns, ISL(:,2), '-*', 'LineWidth', 2);
semilogx(Ns, ISL(:,3), '-d', 'LineWidth', 2);
grid on;
legend('CAN', 'WeCAN', 'G-WeCAN', 'Location', 'northwest');
xlabel('Sequence length N');
ylabel('ISL (dB)');
figure;
semilogx(Ns, time(:,1), '--o', 'LineWidth', 2);
hold on;
semilogx(Ns, time(:,2), '-*', 'LineWidth', 2);
semilogx(Ns, time(:,3), '-d', 'LineWidth', 2);
grid on;
legend('CAN', 'WeCAN', 'G-WeCAN', 'Location', 'northwest');
xlabel('Sequence length N');
ylabel('Elapsed time (sec)');

%%
N = 1000;
M = 3;
alpha = [ones(250,1); zeros(750,1)];
beta = [ones(250,1); zeros(750,1)];
X0 = exp(1i * 2*pi * rand(N, M));
gwecanmimo(N, M, alpha, beta, X0);