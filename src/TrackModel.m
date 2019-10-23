function [dt ErrP ErrF ErrT] = TrackModel(f, K, q, sgn, TMax, f0, dff, dft, WaitBarH, hErrPhi, hErrTau, IFAband, V0, A0, hErrMsg, TrkMode)
% Анализ влияния фильтра на слежение
%    f         - массив дискретизации по частоте
%    K         - комплексный коэффициент передачи
%    sgn       - тип сигнала
%                0 - L1 ПТ
%                1 - L1 ВТ
%                2 - L2 ПТ
%                3 - L2 ВТ
%    N         - количество точек анализа
%    TMax      - диапазон анализа АКФ
%    f0        - центральная частота навигационного сигнала
%    dff       - полоса ФАП
%    dft       - полоса ССЗ
%    WaitBarH  - Waitbar handle
%    hErrPhi   - указатель на график ошибки оценки фазы
%    hErrTau   - указатель на график ошибки оценки задержки
%    IFAband   - полоса УПЧ, Гц (0 - выключен)
%    V0        - начальная скорость, м/с
%    A0        - начальное ускорение, м/с^2
%    hErrMsg   - указатель на индикатор ошибки
%    TrkMode   - режим дальномерных (0) или фазовых (1) измерений
%
% Выходные данные:
%    dt        - шкала дискретизации по времени
%    ErrP      - ошибка по фазе, рад
%    ErrF      - ошибка по частоте, Гц
%    ErrT      - ошибка по задержке, м

global c gd

Gvz = c*mean(gd);


w0 = 2*pi * f0;                     % Несущая частота сигнала

wd = 2*pi * 200e6;                  % Частота дискретизации
Td = 2*pi / wd;                     

wi = wd/4;                          % Промежуточная частота

T = 1e-3;                           % Длительность интевала накопления
Nd = fix(T/Td);
                                    % плотности шума
SgnNd = sqrt(1/4/q/Td);

N = 10;                             % Количество интервалов накопления на каждом этапе

N2 = ceil(TMax/T/N);                % Количество этапов

if (sgn == 0) | (sgn == 2)
    PNSize = 511;   % ГЛОНАСС ПТ
else
    PNSize = 5110;   % ГЛОНАСС ВТ
end

% ПСП
PN = sign(randn(1, PNSize));
TPSP = 1e-3;                        % Длительность ПСП
Tpn = TPSP/PNSize;


% ============================= Параметры ФАП ===============================
SgnA = 1;                         % СКЗ ускорения движения объекта
Sef = 2*0.1*(w0/c)^2*SgnA^2;
SgnEf = sqrt(Sef/T);
SgnNf = sqrt(1/2/q/T*(1+1/2/q/T));
Snf = SgnNf^2*T;
Sdf = (2*q*T)^2;

%dff = 30; % Полоса ФАП
Kf = T * [2*(6/5*dff);  2*(6/5*dff)^2; (6/5*dff)^3];


% ============================ Параметры ССЗ ================================
Set = 2*0.1*SgnA^2;
SgnEt = sqrt(Set/T);
SgnNt = sqrt((c*Tpn)^2*1/4/q/T);
Snt = SgnNt^2*T;
Sdt = 4*q*T/Tpn/c;
% Kt = T * [2*(Set/Snt)^(1/6); 2*(Set/Snt)^(1/3); (Set/Snt)^(1/2)]; % Коэффициетны ССЗ
% dft = 5/6*(Set/Snt)^(1/6);          % Полоса ССЗ
% Set = 1/4*SgnA^2*100;
% SgnEt = sqrt(Set/T);
% SgnNt = sqrt((c*Tpn)^2*1/4/q/T);
% Snt = SgnNt^2*T;
% Sdt = 4*q*T/Tpn/c;
% Kt = T * [ (4*Set/Snt)^0.25; (Set/Snt)^0.5; 0];
% dft = 3/4/sqrt(2)*(Set/Snt)^0.25;

Kt = T * [2*(6/5*dft);  2*(6/5*dft)^2; (6/5*dft)^3];

F = [1 T 0; 0 1 T; 0 0 1];
G = [0; 0; T];

SdQ  =2*q*T;


% Дискретизация по частоте
w = 2*pi*f-w0+wi;
dW = mean(diff(w)); % Шаг моделирования по частоте
NFFT = 2^ceil(log2(2*wd/dW));

% Пересчёт коэффициента передачи в требуемую частотную сетку
wf = (0:NFFT/2-1)/(NFFT/2)*wd/2;
KF = ones(1, NFFT);
KK = conj(K/mean(abs(K)));;
for i=1:size(wf, 2)
    [M I] = min(abs(w-wf(i)));
    if (M <= dW)
        KF(i) = KK(I);
        KF(NFFT-i) = conj(KK(I));
    end
end
KF = conj(KF);

% ======================== Формирование псевдодальности =====================
xf = zeros(3, N*N2);
%xf = [0; 2*pi*100; 2*pi*w0/c*2];
xf = [0; w0/3e8*V0; w0/3e8*A0];
for i=2:N*N2
  xf(:, i) = F*xf(:, i-1) + G*SgnEf*randn(1,1);
end

xt = xf /w0*c;

% %xf0 = zeros(size(xf));
% xf0(:, 1) = xf(:, 1) + 0*[ pi/4; 2*pi*10; 0];
xf0 = xf;

xt0 = xt;
xt0(:, 1) = xt(:, 1) + [-Gvz; 0; 0];

% Udf = zeros(1, N*N2);
% Udt = zeros(1, N*N2);

TStart = tic;

WaitBarH(0, sprintf('%d / %d', 0, N*N2));    

global RcvMdlStopFlag



if IFAband > 0
    if IFAband > wi/2/pi
        IFAband = wi/2/pi;
    end
    [Bi Ai] = ellip(5, 3, 30, (wi + 2*pi*[-IFAband/2 +IFAband/2])/(wd/2));
    Zi = zeros(1, max(length(Ai),length(Bi))-1);
end

for i2=1:N2

  % ======================== Формирование входного процесса ===================
  y = zeros(1, Nd*N);
  for i=1:N
      PNs  = PN(1+mod(fix( ((c+xt(2, (i2-1)*N + i))*Td*(0:Nd-1)+xt(1, (i2-1)*N + i) )/Tpn/c ) , PNSize));
      sgn  = PNs .* exp(1i * ((wi + xf(2, (i2-1)*N + i))*Td*(0:Nd-1) + xf(1, (i2-1)*N + i)));
      obsv = sgn + SgnNd*(randn(1, Nd) + 1i * randn(1, Nd));
      y((i-1)*Nd + (1:Nd)) = obsv;
  end

  if IFAband ~= 0
      % y = filter(Bi, Ai, y);
      [y, Zi] = filter(Bi, Ai, y, Zi);
  end
  
  % Обработка сигнала
  for i=1:fix(size(y, 2)/NFFT)
      y((i-1)*NFFT+(1:NFFT)) = ifft(fft(y((i-1)*NFFT+(1:NFFT))) .* KF);
  end

  

  if i2 == 1
    i0 = 2;
  else
    i0 = 1;
  end
  for i=i0:N
      if RcvMdlStopFlag == 1
          ind = 1:((i2-1)*N+i-1);
          dt = (ind-1)*T;
          ErrP = xf(1, ind) - xf0(1, ind);
          ErrF = (xf(2, ind) - xf0(2, ind))/2/pi;
          ErrT = xt(1, ind) - xt0(1, ind);
          return;
      end
      if (i2-1)*N + i > 1
          done = (i2-1)*N+i;
          WaitBarH(done/N/N2, sprintf('%d / %d  (%8.0f c)', done, N*N2, toc(TStart)/(done-1)*(N*N2-done+1)));
      end
      
    xf0(:, (i2-1)*N + i) = F * xf0(:, (i2-1)*N + i-1);
    %  xf0(:, i) = xf(:, i) + [(i-N/2)/(N/2) * pi; 0; 0];
    xt0(:, (i2-1)*N + i) = F * xt0(:, (i2-1)*N + i-1);
    %    xt0(:, i) = xt(:, i);% + [(i-N/2)/(N/2) * Tpn*c; 0; 0];
    
    PNr  = PN(1+mod(fix( ((c+xt0(2, (i2-1)*N + i))*Td*(0:Nd-1)+xt0(1, (i2-1)*N + i)          )/Tpn/c ) , PNSize));
    PNre = PN(1+mod(fix( ((c+xt0(2, (i2-1)*N + i))*Td*(0:Nd-1)+xt0(1, (i2-1)*N + i) + Tpn/2*c)/Tpn/c ) , PNSize));
    PNrl = PN(1+mod(fix( ((c+xt0(2, (i2-1)*N + i))*Td*(0:Nd-1)+xt0(1, (i2-1)*N + i) - Tpn/2*c)/Tpn/c ) , PNSize));
    
    cr = exp(1i * ((wi + xf0(2, (i2-1)*N + i))*Td*(0:Nd-1) + xf0(1, (i2-1)*N + i)));
   
    I  = 2*SdQ/Nd * real(y((i-1)*Nd + (1:Nd))) * (PNr  .* real(cr))';
    Q  = 2*SdQ/Nd * real(y((i-1)*Nd + (1:Nd))) * (PNr  .* imag(cr))';
    Ie = 2*SdQ/Nd * real(y((i-1)*Nd + (1:Nd))) * (PNre .* real(cr))';
    Il = 2*SdQ/Nd * real(y((i-1)*Nd + (1:Nd))) * (PNrl .* real(cr))';
    
    %    Udf((i2-1)*N + i) = -I*Q / Sdf;
    Udf = -I*Q / Sdf;
    %    Udt((i2-1)*N + i) = sign(I) * (Ie - Il) / Sdt;
    Udt = sign(I) * (Ie - Il) / Sdt;
    
    % xf0(:, (i2-1)*N + i) = xf0(:, (i2-1)*N + i) + Kf * Udf((i2-1)*N + i);  
    % xt0(:, (i2-1)*N + i) = xt0(:, (i2-1)*N + i) + Kt * Udt((i2-1)*N + i);    
    xf0(:, (i2-1)*N + i) = xf0(:, (i2-1)*N + i) + Kf * Udf;  
    xt0(:, (i2-1)*N + i) = xt0(:, (i2-1)*N + i) + Kt * Udt;    
    
  end
  
  
  i = 1:i2*N;
  if (max(abs(xf(1, i)-xf0(1, i))) > pi) | (max(abs(xt(1, i) - xt0(1, i))) > 100)
      set(hErrMsg, 'string', 'Произошёл срыв слежения');
  end
      
  plot(hErrPhi, (i-1)*T, (xf(2, i) - xf0(2, i))/w0*c);
  grid(hErrPhi, 'on');
  
  if TrkMode == 0
      plot(hErrTau, (i-1)*T, xt(1, i) - xt0(1, i));
      grid(hErrTau, 'on');
  else
      plot(hErrTau, (i-1)*T, (xf(1, i) - xf0(1, i))/w0*c + Gvz);
      grid(hErrTau, 'on');
  end

end

dt = (0:N*N2-1)*T;
ErrP = xf(1, :) - xf0(1, :);
ErrF = (xf(2, :) - xf0(2, :))/2/pi;

if TrkMode == 0
    ErrT = xt(1, :) - xt0(1, :);
else
    ErrT = (xf(1, :) - xf0(1, :))/w0*c + Gvz;
end
