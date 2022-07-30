clear all
close all
clc


#número de amostras
n = 100000;

#número de simulações monte carlo
N = 14;

# Valores de Eb/No em db
EbNoDb = [0:N];

#Symbolos por modulação
symbol_qpsk = n/2;
symbol_16QAM = n/4;

#k por modulação
k_16QAM = 4;
k_qpsk = 2;

#gerando vetor aleatório entre 0 e 1
sinal_binario = randi(2, 1, n)-1;

#modulação bpsk
info_bpsk = [];
for i = 1:n
  if sinal_binario(i) == 0;
    info_bpsk(i) = -1;
  else;
    info_bpsk(i) = 1;
  endif
endfor


#modulação BFSK
info_bfsk = [];
for i = 1:n
  if sinal_binario(i) == 0;
    info_bfsk(i) = 1;
  else;
    info_bfsk(i) = j;
  endif
endfor


#modulação QPSK
info_qpsk = [];
sinal_dois_bits = reshape(sinal_binario, 2, [])';
for i = 1:length(sinal_dois_bits(:,1))
  if polyval(sinal_dois_bits(i,:), 2) == 0;
    info_qpsk(i) = 1;
  elseif polyval(sinal_dois_bits(i,:), 2) == 1;
    info_qpsk(i) = j;
  elseif polyval(sinal_dois_bits(i,:), 2) == 2;
    info_qpsk(i) = -j;
  elseif polyval(sinal_dois_bits(i,:), 2) == 3;
    info_qpsk(i) = -1;
  endif
endfor


#modulação 16QAM
info_16QAM = [];
d = sqrt(2/5);
sinal_quatro_bits = reshape(sinal_binario, 4, [])';
for i = 1:length(sinal_quatro_bits(:,1))
  if polyval(sinal_quatro_bits(i,:), 2) == 0;   #0000
    info_16QAM(i) = -1*(3/2)*d + j*(3/2)*d;
  elseif polyval(sinal_quatro_bits(i,:), 2) == 1;   #0001
    info_16QAM(i) = -1*(3/2)*d + j*(1/2)*d;
  elseif polyval(sinal_quatro_bits(i,:), 2) == 2;   #0010
    info_16QAM(i) = -1*(3/2)*d - j*(3/2)*d;
  elseif polyval(sinal_quatro_bits(i,:), 2) == 3;   #0011
    info_16QAM(i) = -1*(3/2)*d - j*(1/2)*d;
  elseif polyval(sinal_quatro_bits(i,:), 2) == 4;   #0100
    info_16QAM(i) = -1*(1/2)*d + j*(3/2)*d;
  elseif polyval(sinal_quatro_bits(i,:), 2) == 5;   #0101
    info_16QAM(i) = -1*(1/2)*d + j*(1/2)*d;
  elseif polyval(sinal_quatro_bits(i,:), 2) == 6;   #0110
    info_16QAM(i) = -1*(1/2)*d - j*(3/2)*d;
  elseif polyval(sinal_quatro_bits(i,:), 2) == 7;   #0111
    info_16QAM(i) = -1*(1/2)*d - j*(1/2)*d;
  elseif polyval(sinal_quatro_bits(i,:), 2) == 8;   #1000
    info_16QAM(i) = 1*(3/2)*d + j*(3/2)*d;
  elseif polyval(sinal_quatro_bits(i,:), 2) == 9;   #1001
    info_16QAM(i) = 1*(3/2)*d + j*(1/2)*d;   
  elseif polyval(sinal_quatro_bits(i,:), 2) == 10;   #1010
    info_16QAM(i) = 1*(3/2)*d - j*(3/2)*d;
  elseif polyval(sinal_quatro_bits(i,:), 2) == 11;   #1011
    info_16QAM(i) = 1*(3/2)*d - j*(1/2)*d;  
  elseif polyval(sinal_quatro_bits(i,:), 2) == 12;   #1100
    info_16QAM(i) = 1*(1/2)*d + j*(3/2)*d;
  elseif polyval(sinal_quatro_bits(i,:), 2) == 13;   #1101
    info_16QAM(i) = 1*(1/2)*d + j*(1/2)*d;
  elseif polyval(sinal_quatro_bits(i,:), 2) == 14;   #1110
    info_16QAM(i) = 1*(1/2)*d - j*(3/2)*d;
  elseif polyval(sinal_quatro_bits(i,:), 2) == 15;   #1111
    info_16QAM(i) = 1*(1/2)*d - j*(1/2)*d;
  endif
endfor


#desvio padrão
EbNo = 10.^(EbNoDb/10);
s_bpsk = sqrt(1./(2*EbNo));
s_bfsk = sqrt(1./(2*EbNo));
s_qpsk = sqrt(1./(2*k_qpsk*EbNo));
s_16QAM = sqrt(1./(2*k_16QAM*EbNo));


#gerando ruídos
ruido_bpsk = [];
ruido_bfsk = [];
ruido_qpsk = [];
ruido_16QAM = [];
for i = 1:length(s_bfsk)
  ruido_bpsk(i,:) = s_bpsk(i)*randn(1, n); #ruído bpsk
  ruido_bfsk(i,:) = s_bfsk(i)*randn(1, n) + s_bfsk(i)*randn(1, n)*j; #ruído bfsk
  ruido_qpsk(i,:) = s_qpsk(i)*randn(1, symbol_qpsk) + s_qpsk(i)*randn(1, symbol_qpsk)*j; #ruído qpsk
  ruido_16QAM(i,:) = s_16QAM(i)*randn(1, symbol_16QAM) + s_16QAM(i)*randn(1, symbol_16QAM)*j; #ruído 16QAM
endfor


#gerando sinal recebido
info_bpsk_rx = ruido_bpsk + info_bpsk;
info_bfsk_rx = ruido_bfsk + info_bfsk;
info_qpsk_rx = ruido_qpsk + info_qpsk;
info_16QAM_rx = ruido_16QAM + info_16QAM;


#Demodulando o sinal bpsk
info_rx_bpsk_demod = [];
for i=1:size(info_bpsk_rx, 1)
  for k=1:size(info_bpsk_rx, 2)
    if(info_bpsk_rx(i, k) <= 0)
      info_rx_bpsk_demod(i, k) = 0;
    else
      info_rx_bpsk_demod(i, k) = 1;
    endif
  endfor
endfor


#Demodulando o sinal bfsk
info_rx_bfsk_demod = [];
for i=1:size(info_bfsk_rx, 1)
  for k=1:size(info_bfsk_rx, 2)
    #convertendo angulo para somente positivo
    angle_symbol = angle(info_bfsk_rx(i, k));
    if(angle_symbol < 0)
      angle_symbol = angle_symbol+(2*pi);
    endif
    
    if(angle_symbol < (pi/4) || angle_symbol > (5/4)*pi)
      info_rx_bfsk_demod(i, k) = 0;
    else
      info_rx_bfsk_demod(i, k) = 1;
    endif
  endfor
endfor


#Demodulando o sinal QPSK
info_rx_qpsk_demod = [];
info_rx_qpsk = [];
for i=1:size(info_qpsk_rx, 1)
  qpsk_demod = [];
  for k=1:size(info_qpsk_rx, 2)
    #convertendo angulo para somente positivo
    angle_symbol = angle(info_qpsk_rx(i, k));
    if(angle_symbol < 0)
      angle_symbol = angle_symbol+(2*pi);
    endif

    if(angle_symbol <= (1/4)*pi || angle_symbol > (7/4)*pi)
      qpsk_demod = [qpsk_demod, 0, 0];
      info_rx_qpsk(i, k) = 1;
    elseif(angle_symbol <= (3/4)*pi && angle_symbol > (1/4)*pi)
      qpsk_demod = [qpsk_demod, 0, 1];
      info_rx_qpsk(i, k) = j;
    elseif(angle_symbol <= (5/4)*pi && angle_symbol > (3/4)*pi)
      qpsk_demod = [qpsk_demod, 1, 1];
      info_rx_qpsk(i, k) = -1;
    elseif(angle_symbol <= (7/4)*pi && angle_symbol > (5/4)*pi)
      qpsk_demod = [qpsk_demod, 1, 0];
      info_rx_qpsk(i, k) = -j;
    endif
  endfor
  info_rx_qpsk_demod = [info_rx_qpsk_demod; qpsk_demod];
endfor


#Demodulando o sinal 16QAM
info_rx_16AQAM_demod = [];
info_rx_16AQAM = [];
for i=1:size(info_16QAM_rx, 1)
  qam_demod = [];
  for k=1:size(info_16QAM_rx, 2)
    if(real(info_16QAM_rx(i, k)) > d)
      if(imag(info_16QAM_rx(i, k)) > d)
        qam_demod = [qam_demod, 1, 0, 0, 0];
        info_rx_16AQAM(i,k) = 1*(3/2)*d + j*(3/2)*d;
      elseif(imag(info_16QAM_rx(i, k)) > 0) && (imag(info_16QAM_rx(i, k)) <= d)
        qam_demod = [qam_demod, 1, 0, 0, 1];
        info_rx_16AQAM(i,k) = 1*(3/2)*d + j*(1/2)*d; 
      elseif(imag(info_16QAM_rx(i, k)) > -d) && (imag(info_16QAM_rx(i, k)) <= 0)
        qam_demod = [qam_demod, 1, 0, 1, 1];
        info_rx_16AQAM(i,k) = 1*(3/2)*d - j*(1/2)*d; 
      elseif(imag(info_16QAM_rx(i, k)) <= -d)
        qam_demod = [qam_demod, 1, 0, 1, 0];
        info_rx_16AQAM(i,k) = 1*(3/2)*d - j*(3/2)*d;
      endif
    elseif(real(info_16QAM_rx(i, k)) >=0) && (real(info_16QAM_rx(i, k)) < d)
      if(imag(info_16QAM_rx(i, k)) > d)
        qam_demod = [qam_demod, 1, 1, 0, 0];
        info_rx_16AQAM(i,k) = 1*(1/2)*d + j*(3/2)*d;
      elseif(imag(info_16QAM_rx(i, k)) > 0) && (imag(info_16QAM_rx(i, k)) <= d)
        qam_demod = [qam_demod, 1, 1, 0, 1];
        info_rx_16AQAM(i,k) = 1*(1/2)*d + j*(1/2)*d;
      elseif(imag(info_16QAM_rx(i, k)) > -d) && (imag(info_16QAM_rx(i, k)) <= 0)
        qam_demod = [qam_demod, 1, 1, 1, 1];
        info_rx_16AQAM(i,k) = 1*(1/2)*d - j*(1/2)*d;
      elseif(imag(info_16QAM_rx(i, k)) <= -d)
        qam_demod = [qam_demod, 1, 1, 1, 0];
        info_rx_16AQAM(i,k) = 1*(1/2)*d - j*(3/2)*d;
      endif
    elseif(real(info_16QAM_rx(i, k)) >= -d) && (real(info_16QAM_rx(i, k)) < 0)
      if(imag(info_16QAM_rx(i, k)) > d)
        qam_demod = [qam_demod, 0, 1, 0, 0];
        info_rx_16AQAM(i,k) = -1*(1/2)*d + j*(3/2)*d;
      elseif(imag(info_16QAM_rx(i, k)) > 0) && (imag(info_16QAM_rx(i, k)) <= d)
        qam_demod = [qam_demod, 0, 1, 0, 1];
        info_rx_16AQAM(i,k) = -1*(1/2)*d + j*(1/2)*d;
      elseif(imag(info_16QAM_rx(i, k)) > -d) && (imag(info_16QAM_rx(i, k)) <= 0)
        qam_demod = [qam_demod, 0, 1, 1, 1];
        info_rx_16AQAM(i,k) = -1*(1/2)*d - j*(1/2)*d;
      elseif(imag(info_16QAM_rx(i, k)) <= -d)
        qam_demod = [qam_demod, 0, 1, 1, 0];
        info_rx_16AQAM(i,k) = -1*(1/2)*d - j*(3/2)*d;
      endif
    elseif(real(info_16QAM_rx(i, k)) < -d)
      if(imag(info_16QAM_rx(i, k)) > d)
        qam_demod = [qam_demod, 0, 0, 0, 0];
        info_rx_16AQAM(i,k) = -1*(3/2)*d + j*(3/2)*d;
      elseif(imag(info_16QAM_rx(i, k)) > 0) && (imag(info_16QAM_rx(i, k)) <= d)
        qam_demod = [qam_demod, 0, 0, 0, 1];
        info_rx_16AQAM(i,k) = -1*(3/2)*d + j*(1/2)*d;
      elseif(imag(info_16QAM_rx(i, k)) > -d) && (imag(info_16QAM_rx(i, k)) <= 0)
        qam_demod = [qam_demod, 0, 0, 1, 1];
        info_rx_16AQAM(i,k) = -1*(3/2)*d - j*(1/2)*d;
      elseif(imag(info_16QAM_rx(i, k)) <= -d)
        qam_demod = [qam_demod, 0, 0, 1, 0];
        info_rx_16AQAM(i,k) = -1*(3/2)*d - j*(3/2)*d;
      endif
    endif
  endfor
  info_rx_16AQAM_demod = [info_rx_16AQAM_demod; qam_demod];
endfor

  
#cálculo erros para o sinal bpsk
erros_bpsk = (info_rx_bpsk_demod != sinal_binario);
erros_bfsk = (info_rx_bfsk_demod != sinal_binario);
erros_qpsk = (info_rx_qpsk_demod != sinal_binario);
erros_16QAM = (info_rx_16AQAM_demod != sinal_binario);
pb_erro_bpsk = [];
pb_erro_bfsk = [];
pb_erro_qpsk = [];
pb_erro_16QAM = [];
for i=1:size(erros_bpsk, 1)
  pb_erro_bpsk(i) = sum(erros_bpsk(i,:) == 1)/n;
  pb_erro_bfsk(i) = sum(erros_bfsk(i,:) == 1)/n;
  pb_erro_qpsk(i) = sum(erros_qpsk(i,:) == 1)/n;
  pb_erro_16QAM(i) = sum(erros_16QAM(i,:) == 1)/n;
endfor

#cálculo erro symbol qpsk
pb_erro_symbol_qpsk = [];
erros_symbol_qpsk = info_rx_qpsk != info_qpsk;
for i=1:size(erros_symbol_qpsk, 1)
    pb_erro_symbol_qpsk(i) = sum(erros_symbol_qpsk(i,:) == 1)/(n/2);
endfor

#cálculo erro symbol 16QAM
pb_erro_symbol_16QAM = [];
erros_symbol_16QAM = info_rx_16AQAM != info_16QAM;
for i=1:size(erros_symbol_16QAM, 1)
    pb_erro_symbol_16QAM(i) = sum(erros_symbol_16QAM(i,:) == 1)/(n/4);
endfor

#cálculos teóricos
#erros de bit
pb_erro_bfsk_teor = qfunc(sqrt(EbNo));
pb_erro_bpsk_teor = qfunc(sqrt(2*EbNo));
pb_erro_qpsk_teor = (2/log2(4))*(qfunc (sqrt(2*EbNo*log2(4))*sin(pi/4)));
pb_erro_16QAM_teor = (3/4)*(qfunc (sqrt(EbNo*4/5)));
#erros de símbolos
pb_erro_qpsk_symbol_teor = 2*(qfunc (sqrt(4*EbNo)*sin(pi/4)));
pb_erro_16QAM_symbol_teor = 3*(qfunc (sqrt(4*EbNo*1/5)));

#plotando gráfico BER
figure(1)
semilogy(EbNoDb, pb_erro_bfsk, "g")
hold on
semilogy(EbNoDb, pb_erro_bpsk, "m")
semilogy(EbNoDb, pb_erro_qpsk, "c")
semilogy(EbNoDb, pb_erro_16QAM, "b")
semilogy(EbNoDb, pb_erro_bfsk_teor, "**g")
semilogy(EbNoDb, pb_erro_bpsk_teor, "++m")
semilogy(EbNoDb, pb_erro_qpsk_teor, "xc")
semilogy(EbNoDb, pb_erro_16QAM_teor, "--b")
t = title ("Probabilidade de Erro de Bit - BER");
set (t, "fontsize", 14);
h = legend ("BFSK", "BPSK", "QPSK", "16QAM", "BFSK - Teor", "BPSK - Teor", "QPSK-Teor", "16AQAM - Teor" );
legend (h, "location", "northeastoutside");
set (h, "fontsize", 12);
grid on;
ylabel("Probabilidade de erro");
xlabel("EB/No");
xlim([0 15])
ylim([10^-5 1])

#plotando gráfico SER
figure(2)
semilogy(EbNoDb, pb_erro_symbol_qpsk, "c")
hold on
grid on;
semilogy(EbNoDb, pb_erro_symbol_16QAM, "b")
semilogy(EbNoDb, pb_erro_qpsk_symbol_teor, "--c")
semilogy(EbNoDb, pb_erro_16QAM_symbol_teor, "--b")
t = title ("Probabilidade de Erro de Símbolo - SER");
set (t, "fontsize", 14);
h = legend ("QPSK", "16QAM", "QPSK - Teor", "16QAM - Teor");
legend (h, "location", "northeastoutside");
set (h, "fontsize", 12);
grid on;
ylabel("Probabilidade de erro");
xlabel("EB/No");
xlim([0 15])
ylim([10^-5 1])


