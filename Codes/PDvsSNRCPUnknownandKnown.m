%PD vs SNR CP unknown and known
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
SNR_dB = linspace(-25,-0,Td);
Pfa = 0.05;
rho1 = zeros(length(SNR_dB),1);
for i=1:length(SNR_dB)
    rho1(i) = (Tc)/(Td+Tc)*(db_to_lin(SNR_dB(i)))/(1+db_to_lin(SNR_dB(i)));
end
 etal =(1/sqrt(M)).*erfcinv(2.*Pfa);
 Pdcpunknowntheory = zeros(length(SNR_dB),1);
 for i=1:length(SNR_dB)
     Pdcpunknowntheory(i) = 0.5.*erfc(sqrt(M).*(etal - rho1(i)/(1-rho1(i) .^2)));
 end
 %%
rhoc = zeros(length(SNR_dB),1);
for i=1:length(SNR_dB)
    rhoc(i) = (db_to_lin(SNR_dB(i)))/(1+db_to_lin(SNR_dB(i)));
end
 etal =(1/sqrt(M1)).*erfcinv(2.*Pfa);
 Pdcpknowntheory = zeros(length(SNR_dB),1);
 for i=1:length(SNR_dB)
     Pdcpknowntheory(i) = 0.5.*erfc(sqrt(M1).*(etal - rhoc(i)/(1-rhoc(i) .^2)));
 end
 %%
figure
plot(SNR_dB,Pdcpknowntheory,SNR_dB,Pdcpunknowntheory)
grid on
xlabel(" SNR (dB)",'interpreter','latex');
ylabel("Probability of Detection ($P_D$)",'interpreter','latex');
title(" $P_d$ vs SNR for $ P_{fa} = -10dB$ CP known and unknown ",'interpreter','latex');
legend('CP Known','CP unknown')