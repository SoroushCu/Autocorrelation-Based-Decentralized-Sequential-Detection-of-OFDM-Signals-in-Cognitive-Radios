%ROC plot for fixed SNR AWGN and Theory
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
SNR_dB = -10;
SNR = db_to_lin(SNR_dB);
rho1 = (Tc)/(Td+Tc)*((SNR)/(SNR+1));
Pfa = linspace(-1,1,1000);
etal = zeros(length(Pfa),1);
for i=1:length(Pfa)
    etal(i) =(1/sqrt(M)).*erfcinv(2.*Pfa(i));
end
Pd = 0.5.*erfc(sqrt(M).*(etal - rho1/(1-rho1 .^2)));

figure
plot(Pfa,Pd)
grid on
xlabel(" Probability of False Alarm ($P_{fa}$)",'interpreter','latex');
ylabel("Probability of Detection ($P_D$)",'interpreter','latex');
title(" ROC Plot for $ SNR = -10dB$ ",'interpreter','latex');
%%
Nc = 100;       
Td = 32;          
Tc = 8;     
Ts = Td + Tc;
N = Td * Nc;
Mod = 16;                 
k = log2(Mod);
SNR_dB = -10;        
Pfasim = linspace(0,1,1000);
Pdsim = zeros(1000);
data = randi([0, 1], k*N, 1);
dataSymbols = reshape(data, k, []).';
constellation = qammod(0:Mod-1, Mod);
x = constellation(bin2dec(num2str(dataSymbols))+1);
x = reshape(x, Td, Nc);
xifft = ifft(x, Td);
xifftwithcp = [xifft(end-Tc+1:end, :); xifft];
xsend = xifftwithcp(:);
N0 = 1/(2*db_to_lin(SNR_dB));
thresh = erfcinv(2.*Pfasim).*1.3*N0+2.9;
for i = 1:length(Pfasim)
    noise = sqrt(N0/2) * (randn(size(xsend)) + 1j * randn(size(xsend)));
    y = xsend + noise;
    yCP = reshape(y, Ts, []);
    y = yCP(Tc+1:end, :);
    Y = fft(y, Td);
    Y = Y(:);
    Det = sum(abs(Y) > thresh(i));
    Pdsim(i) = Det / (N);
end
disp('shit')
figure
plot(Pfa,Pd)
hold on
plot(Pfasim, Pdsim);
xlabel(" Probability of False Alarm ($P_{fa}$)",'interpreter','latex');
ylabel("Probability of Detection ($P_D$)",'interpreter','latex');
title(" ROC Plot for $ SNR = -10dB$ ",'interpreter','latex');
legend('Theory','Simulation')
grid on;