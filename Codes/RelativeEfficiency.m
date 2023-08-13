%Relative Efficiency
%%
clc
clear all
Td = 32;
Tc = 8;
Ts = Tc+Td;
M = 100*(Td+Tc);
M1 = Tc*M/Ts;
OFDMBW = 5e6;
Nc = 100;%Number of OFDM blocks
%%
SNR_dB = linspace(-20,-10,Td);
Pd = 0.95;
Pfa = 1-Pd;
beta = Pfa;
A = (1-beta)/Pfa;
B = beta/(1-Pfa);
rho1 = zeros(length(SNR_dB),1);
for i=1:length(SNR_dB)
    rho1(i) = (Tc)/(Td+Tc)*(db_to_lin(SNR_dB(i)))/(1+db_to_lin(SNR_dB(i)));
end
ELn_H0 =  zeros(length(SNR_dB),1);
ELn_H1 =  zeros(length(SNR_dB),1);
for i=1:length(SNR_dB)
    ELn_H0(i) = -M.*log10(1-rho1(i)^2)-2*M*((rho1(i)^2)/(1-rho1(i) ^ 2));
    ELn_H1(i) = -M.*log10(1-rho1(i)^2);
end

EKs_H0 =  zeros(length(SNR_dB),1);
EKs_H1 =  zeros(length(SNR_dB),1);

for i=1:length(SNR_dB)
   EKs_H0(i) = (Pfa.*log10(A)+(1-Pfa).*log10(B))./( ELn_H0(i));
   EKs_H1(i) = ( (1-beta).*log10(A) + beta.*log10(B))/(ELn_H1(i));
end
Km = max( EKs_H0,EKs_H1);

rhon = rho1;
Kf = zeros(length(SNR_dB),1);
for i=1:length(SNR_dB)
    Kf(i) = (erfcinv(2*Pfa) - (1-rho1(i).^2).*erfcinv(2*Pd) ).^2./(M*rho1(i) .^ 2);
end
RE = Kf./Km(:,1);

figure
plot(SNR_dB,RE,'r--')
grid on
xlabel(" SNR (dB)",'interpreter','latex');
ylabel("Average Sample Number (ASN)",'interpreter','latex');
title(" $P_{fa} = 0.05$ , $P_D = 0.95$ ",'interpreter','latex');
ylim([0,4])