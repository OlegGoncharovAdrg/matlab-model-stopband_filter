function [t s y PSP] = RfAffectAnalyse(f, K, sgn, N, TauMax, f0, WaitBarH, IFAband)
%path(path, '../../../Sub');
% Анализ искажения автокорреляционной функции
%    f         - массив дискретизации по частоте
%    K         - комплексный коэффициент передачи
%    sgn       - тип сигнала
%                0 - L1 ПТ
%                1 - L1 ВТ
%                2 - L2 ПТ
%                3 - L2 ВТ
%    N         - количество точек анализа
%    TauMax    - диапазон анализа АКФ
%    f0        - центральная частота навигационного сигнала
%    WaitBarH  - Waitbar handle
%    IFAband   - полоса фильтра УПЧ, Гц (0 - отключить)
%
% Выходные данные:
%    t     - шкала дискретизации по времени
%    s     - чистый сигнал
%    y     - искажённые сигнал
    

%w0 = 2*pi*f0;
w0 = 2*pi * 1602e6;

wd = 2*pi*4000e6;
Td = 2*pi/wd;

wi = wd/4;

randn('seed', 0);

% BOC
if (sgn == 0) | (sgn == 2)
    PNSize = 1023;
    PN = sign(randn(1, PNSize));
    PNSize = 2*PNSize;
    PN = expand_arr1(PN, PNSize) .* (1-2*(mod(0:PNSize-1, 2)));
else
    PNSize = 5110;
    PN = sign(randn(1, PNSize));
    PNSize = 2*PNSize;
    PN = expand_arr1(PN, PNSize) .* (1-2*(mod(0:PNSize-1, 2)));
end
Tpn = 1e-3/PNSize;
Tpsp = Tpn;



% % BOC
% m = 1;
% n = 7;
% PN = sign(randn(1, fix(m*1023)));
% M = -1+2*mod(1:fix(n*1023), 2);
% PNSize = max(m, 2*n)*1023;
% Tpn = 1e-3/PNSize;
% PN = expand_arr1(PN, PNSize) .* expand_arr1(M, PNSize);
% Tpsp = 1e-3/m/1023;


%T = 1e-3;
T = TauMax;
Nd = fix(T/Td);
%Nd = N;

% Дискретизация по времени
dt = ((0:N-1)-N/2+0.5)/(N/2) * TauMax;


% Перенос сигнала на центральную частоту анализа
w = 2*pi*f-w0+wi;

% Параметры БПФ
dW = mean(diff(w));
NFFT = 2^ceil(log2(2*wd/dW));

% Формирование анализируемого сигнала
ref = exp(1i*wi*Td*(0:Nd+NFFT-1));
PSP = PN(1+mod(fix((Td*(0:Nd+NFFT-1))/Tpn), PNSize));
s = ref .* PSP;
y = s;

% Частоты анализа БПФ
wf = (0:NFFT/2-1)/(NFFT/2)*wd/2;

% Распределение отсчётов ЧХ по соответствующим каналам анализа БПФ
KF = ones(1, NFFT);

% Нормировка коэффициента передачи в пределах полосы сигнала
ind = find((f>f0-1/Tpn) & (f<f0+1/Tpn));
KK = conj(K/mean(abs(K(ind))));  

for i=1:size(wf, 2)
    [M I] = min(abs(w-wf(i)));
    if (M <= dW)
        KF(i) = KK(I);
        KF(NFFT-i) = conj(KK(I));
    end
end
KF = conj(KF);

if IFAband > 0
    if IFAband > wi/2/pi
        IFAband = wi/2/pi;
    end
    [Bi Ai] = ellip(5, 3, 30, (wi + 2*pi*[-IFAband/2 +IFAband/2])/(wd/2));
    Zi = zeros(1, max(length(Ai), length(Bi))-1);
    y1 = filter(Bi, Ai, y);
    ref1 = filter(Bi, Ai, ref);
else
    y1 = y;
    ref1 = ref;
end


% Искажённый в фильтре сигнал
y2 = zeros(size(y));
for i=1:fix(size(y, 2)/NFFT)
    y2((i-1)*NFFT+(1:NFFT)) = ifft(fft(y1((i-1)*NFFT+(1:NFFT))) .* KF);
    ref1((i-1)*NFFT+(1:NFFT)) = ifft(fft(ref((i-1)*NFFT+(1:NFFT))) .* KF);    
end

y = y(1:Nd);       % Сигнал не искажённый
y1 = y1(1:Nd);     % Сигнал искажённый в УПЧ
y2 = y2(1:Nd);

s = s(1:Nd);
y = y2(1:Nd);
t = Td*(0:Nd-1);
PSP = PSP(1:Nd);

% ref = ref(1:Nd);   % Сигнал опорный (гармонический)
% ref1 = ref1(1:Nd);   % Сигнал опорный (гармонический)

% % Расчёт корреляционной функции
% R = zeros(1, N);
% R1 = zeros(1, N);


% TStart = tic;
% for i=1:N
%     if (i>1)
%         %        waitbar(i/N, WaitBarH, sprintf('%3d%% (%5d c)', ceil(100*i/N), ceil(toc(TStart)/(i-1)*(N-i+1))));
%         WaitBarH(i/N, sprintf('%3d%% (%5d c)', ceil(100*i/N), ceil(toc(TStart)/(i-1)*(N-i+1))));
%     end
    
%     PNr = PN(1+mod(fix((-dt(i)+Td*(0:Nd-1))/Tpn), PNSize));
%     %          PN(1+fix(mod((      Td*(0:Nd+NFFT-1))/Tpn, PNSize)))
    
%     R(i)  = abs(y  * (ref  .* PNr)');
%     R1(i) = abs(y2 * (ref1 .* PNr)');
%     % R(i)  = abs(y  * (ref .* PNr)');
%     % R1(i) = abs(y2 * (ref .* PNr)');
% end

% MX = max(R);
% % R = R(N:-1:1)/MX;
% % R1 = R1(N:-1:1)/MX;

