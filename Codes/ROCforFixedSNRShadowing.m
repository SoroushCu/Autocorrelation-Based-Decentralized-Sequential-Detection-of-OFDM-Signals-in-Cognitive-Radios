%ROC plot for fixed SNR with Shadowing
%%
clc
clear all
Td = 32;
Tc = 8;
Ts = Tc+Td;
M = 100*(Td+Tc);
M1 = Tc*M/Ts;
OFDMBW = 5e6;
Nc = 100;
%%
Nc = 100;       
Td = 32;          
Tc = 8;     
Ts = Td + Tc;
N = Td * Nc;
Mod = 16;                 
k = log2(Mod);
SNR_dB = -10;        
Pfasim2 = linspace(0,1,1000);
Pdsim2 = zeros(1000);
data = randi([0, 1], k*N, 1);
dataSymbols = reshape(data, k, []).';
constellation = qammod(0:Mod-1, Mod);
x = constellation(bin2dec(num2str(dataSymbols))+1);
x = reshape(x, Td, Nc);
xifft = ifft(x, Td);
xifftwithcp = [xifft(end-Tc+1:end, :); xifft];
xsend = xifftwithcp(:);
Shadowstd = 5;
ShadowdB = randn(size(xsend)) * Shadowstd;
Shadow = 10.^(ShadowdB/10);
N0 = 1/(2*db_to_lin(SNR_dB));
thresh = erfcinv(2.*Pfasim2).*1.3*N0+2.9;
for i = 1:length(Pfasim2)
    noise = sqrt(N0/2) * (randn(size(xsend)) + 1j * randn(size(xsend)));
    y = xsend.*Shadow + noise;
    yCP = reshape(y, Ts, []);
    y = yCP(Tc+1:end, :);
    Y = fft(y, Td);
    Y = Y(:);
    Det = sum(abs(Y) > thresh(i));
    Pdsim2(i) = Det / (N);
end
figure;
plot(Pfasim2,Pdsim2)
xlabel(" Probability of False Alarm ($P_{fa}$)",'interpreter','latex');
ylabel("Probability of Detection ($P_D$)",'interpreter','latex');
title(" ROC Plot for $ SNR = -10dB$ shadowing included ",'interpreter','latex');
grid on;