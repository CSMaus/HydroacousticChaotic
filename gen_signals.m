%%%%% Генерация сигналов для разных систем и их исследование
clear all
xx1 = 30;
xx2 = 7;
flag = 1;
[Fsamp,r,FLAG]=mydlg(xx1,xx2,flag); %Загрузка окны выбора частоты и расстояния

%%%%%%%%%%%%%%%%%%%%%% Proto-Lorenz %%%%%%%%%%%%%%%%%%%%%%%%%
tendL = 70; %Время
hL = 0.001; %шаг

%НУ
L_x0 = 1.287;
L_y0 = -4.128;
L_z0 = 41.42;

%Параметры системы
La=10;
Lb=8/3;
Lc=28;

fL = @(t,Lx) [1/(2*(Lx(1)^2+Lx(2)^2))*(-La*Lx(1)^3 + (2*La+Lc-Lx(3))*Lx(1)^2*Lx(2)+(La-2)*Lx(1)*Lx(2)^2-(Lc-Lx(3))*Lx(2)^3);...
    1/(2*(Lx(1)^2+Lx(2)^2))*((Lc-Lx(3))*Lx(1)^3+(La-2)*Lx(1)^2*Lx(2)+(-2*La-Lc+Lx(3))*Lx(1)*Lx(2)^2-La*Lx(2)^3);...
    2*Lx(1)^3*Lx(2)-2*Lx(1)*Lx(2)^3-Lb*Lx(3)];
[t,Lx] = ode45(fL,[0:hL:tendL],[L_x0 L_y0 L_z0]);

Lx(:,1) = Lx(:,1)./max(abs(Lx(:,1)));
Lx(:,2) = Lx(:,2)./max(abs(Lx(:,2)));
Lx(:,3) = Lx(:,3)./max(abs(Lx(:,3)));
    figure('Name','Proto-Lorenz')
    plot3(Lx(:,1),Lx(:,2),Lx(:,3)) %Построение графика 3d
    grid on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
%     set(gca, 'XTick', [-0.5 0 0.5], 'YTick', [-0.5 0 0.5],'ZTick', [0 0.25 0.5 0.75])
PLorenz = Lx(:,1)./max(abs(Lx(:,1))); %Нормирую вектор x к максимуму
PLorenz = transpose(PLorenz);
n = 1:length(PLorenz);
figure('Name','Signals PLorenz')
plot(n/1000, PLorenz)
xlabel('[n] \times 1000, отсчетов')
ylabel('Нормированная амплитуда')
set(gca,'FontName','Times New Roman','FontSize',20)
grid on


%% %%%%%%%%%%%%%%%%%%%%%% Sprott Case B %%%%%%%%%%%%%%%%%%%%%%
% waitbar(0.14,f,'Calculating Sprott Case B...');
tendSp = tendL*4; %частота у Спротта больше чем у Лоренца в 5 раз (по сравнению с изначальной)
hSp = hL*4;
% tendSp = tendL;
% hSp = hL;

% %%%%%%%%%% Case A %%%%%%%%%%%

%НУ
% Sp_Ax0 = 1.3379;
% Sp_Ay0 = 0.9376;
% Sp_Az0 = 0.8454;
Sp_Ax0 = 0.25;
Sp_Ay0 = 0.2625;
Sp_Az0 = 0.18;

tendSpA = tendL*6;
hSpA = hL*6;

fS = @(t,SpAx) [SpAx(2); -SpAx(1) + SpAx(2)*SpAx(3); 1 - SpAx(2)*SpAx(2)];
[t,SpAx] = ode45(fS,[0:hSpA:tendSpA],[Sp_Ax0 Sp_Ay0 Sp_Az0]);
SpAx(:,1) = SpAx(:,1)./max(abs(SpAx(:,1)));
SpAx(:,2) = SpAx(:,2)./max(abs(SpAx(:,2)));
SpAx(:,3) = SpAx(:,3)./max(abs(SpAx(:,3)));
%     figure(3)
%     plot3(SpAx(:,1), SpAx(:,2), SpAx(:,3))
%     title('Sprott Case A')
%     grid on
SprottA = SpAx(:,1)./max(abs(SpAx(:,1))); %Нормирую вектор x к максимуму
SprottA = transpose(SprottA);
figure('Name','Signals SprottA')
plot(n/1000, SprottA)
xlabel('[n] \times 1000, отсчетов')
ylabel('Нормированная амплитуда')
set(gca,'FontName','Times New Roman','FontSize',20)
grid on

%%%%%%%%%% Case B %%%%%%%%%%%
% Sp_Bx0 = 0.7290;
% Sp_By0 = 0.9814;
% Sp_Bz0 = 0.1978;
Sp_Bx0 = 0.4079;
Sp_By0 = 0.7376;
Sp_Bz0 = 0.5345;

fS = @(t,SpBx) [SpBx(2)*SpBx(3); SpBx(1) - SpBx(2); 1 - SpBx(1)*SpBx(2)];
[t,SpBx] = ode45(fS,[0:hSp:tendSp],[Sp_Bx0 Sp_By0 Sp_Bz0]);
SpBx(:,1) = SpBx(:,1)./max(abs(SpBx(:,1)));
SpBx(:,2) = SpBx(:,2)./max(abs(SpBx(:,2)));
SpBx(:,3) = SpBx(:,3)./max(abs(SpBx(:,3)));
%     figure(2)
%     plot3(SpBx(:,1), SpBx(:,2), SpBx(:,3))
%     title('Sprott Case B')
%     grid on
SprottB = SpBx(:,1)./max(abs(SpBx(:,1))); %Нормирую вектор x к максимуму
SprottB = transpose(SprottB);
figure('Name','Signals SprottB')
plot(n/1000, SprottB)
xlabel('[n] \times 1000, отсчетов')
ylabel('Нормированная амплитуда')
set(gca,'FontName','Times New Roman','FontSize',20)
grid on

%%%%%%%%%% Case C %%%%%%%%%%%
Sp_Cx0 = 0.4079;
Sp_Cy0 = 0.7376;
Sp_Cz0 = 0.5454;
fS = @(t,SpCx) [SpCx(2)*SpCx(3); SpCx(1) - SpCx(2); 1 - SpCx(1)*SpCx(1)];
[t,SpCx] = ode45(fS,[0:hSp:tendSp],[Sp_Cx0 Sp_Cy0 Sp_Cz0]);
SpCx(:,1) = SpCx(:,1)./max(abs(SpCx(:,1)));
SpCx(:,2) = SpCx(:,2)./max(abs(SpCx(:,2)));
SpCx(:,3) = SpCx(:,3)./max(abs(SpCx(:,3)));
%     figure(4)
%     plot3(SpCx(:,1), SpCx(:,2), SpCx(:,3))
%     title('Sprott Case C')
%     grid on
SprottC = SpCx(:,1)./max(abs(SpCx(:,1))); %Нормирую вектор x к максимуму
SprottC = transpose(SprottC);
figure('Name','Signals SprottC')
plot(n/1000, SprottC)
xlabel('[n] \times 1000, отсчетов')
ylabel('Нормированная амплитуда')
set(gca,'FontName','Times New Roman','FontSize',20)
grid on

%%%%%%%%%%% Case E %%%%%%%%%%%
% Sp_Ex0 = 0.5491;
% Sp_Ey0 = 0.3614;
% Sp_Ez0 = 0.1878;
Sp_Ex0 = 0.215;
Sp_Ey0 = 0.3625;
Sp_Ez0 = 0.96;
fS = @(t,SpEx) [SpEx(2)*SpEx(3); SpEx(1)*SpEx(1) - SpEx(2); 1 - 4*SpEx(1)];
[t,SpEx] = ode45(fS,[0:hSp:tendSp],[Sp_Ex0 Sp_Ey0 Sp_Ez0]);
SpEx(:,1) = SpEx(:,1)./max(abs(SpEx(:,1)));
SpEx(:,2) = SpEx(:,2)./max(abs(SpEx(:,2)));
SpEx(:,3) = SpEx(:,3)./max(abs(SpEx(:,3)));
%     figure(5)
%     plot3(SpEx(:,1), SpEx(:,2), SpEx(:,3))
%     title('Sprott Case E')
%     xlabel('X')
%     ylabel('y')
%     zlabel('z')
%     grid on
SprottE = SpEx(:,1)./max(abs(SpEx(:,1))); %Нормирую вектор x к максимуму
SprottE = transpose(SprottE);
figure('Name','Signals SprottE')
plot(n/1000, SprottE)
xlabel('[n] \times 1000, отсчетов')
ylabel('Нормированная амплитуда')
set(gca,'FontName','Times New Roman','FontSize',20)
grid on


figure('Name','Sprott systems')
subplot(2,2,1)
    plot3(SpAx(:,1), SpAx(:,2), SpAx(:,3))
    title('Sprott Case A')
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    grid on
subplot(2,2,2)
    plot3(SpBx(:,1), SpBx(:,2), SpBx(:,3))
    title('Sprott Case B')
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    grid on
subplot(2,2,3)
    plot3(SpCx(:,1), SpCx(:,2), SpCx(:,3))
    title('Sprott Case C')
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    grid on
subplot(2,2,4)
    plot3(SpEx(:,1), SpEx(:,2), SpEx(:,3))
    title('Sprott Case E')
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    grid on


%% %%%%%%%%%%%% Four-wing (Fw) chaotic attractor new 3D system %%%%%%%%%
%%%A four-wing chaotic attractor generated from a new
%%%3-D quadratic autonomous system Guoyuan Qi, Guanrong Chen etc.
% waitbar(0.28,f,'Calculating Four-wing new 3D sys...');

tendFw = tendL/4;
hFw = hL/4;
% tendFw = tendL;
% hFw = hL;


%Params

% a = 1.6;
% b = 3.34;
% c = -0.9;
% d = 1.5;
% e = 0.25;
a = 19;
b = 35;
c = -1;
d = 16;
e = 3;

f = c*e -a;
g = sqrt((a+c*e)+4*a*d*e);
h = f - g - 2*c*e;
i = f + g - 2*c*e;

%НУ
xFw = sqrt(b*d*(f-g)/h);
zFw = (f - g)/(2*e);
yFw = (zFw - c)*xFw/d;

fFw = @(t,Fwx) [a*(Fwx(2) - Fwx(1)) + e*Fwx(2)*Fwx(3); c*Fwx(1) + d*Fwx(2) - Fwx(1)*Fwx(3); -b*Fwx(3) + Fwx(1)*Fwx(2)];
[t,Fwx] = ode45(fFw,[0:hFw:tendFw],[xFw yFw zFw]);


Fwx(:,1) = Fwx(:,1)./max(abs(Fwx(:,1)));
Fwx(:,2) = Fwx(:,2)./max(abs(Fwx(:,2)));
Fwx(:,3) = Fwx(:,3)./max(abs(Fwx(:,3)));
    figure('Name','Four wing by new 3D system')
    plot3(Fwx(:,1), Fwx(:,2), Fwx(:,3))
    grid on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
Four_wing = Fwx(:,1)./max(abs(Fwx(:,1)));
Four_wing = transpose(Four_wing);
figure('Name','Signals Four_wing')
plot(n/1000, Four_wing)
xlabel('[n] \times 1000, отсчетов')
ylabel('Нормированная амплитуда')
set(gca,'FontName','Times New Roman','FontSize',20)
grid on

%%%%%%%%%%%%% Знакопеременный чирп %%%%%%%%%%%%
% waitbar(0.42,f,'Calculating Chirp...');
for i = 1:tendL/hL+1
    tspan = i*hL;
    w = 1.2;%*(tspan/100)для дополнтельной частотной модуляции чирпа
    Amax = 1.7; %масимальная амплитуда синуса, которую он достигает к конце периода
    Tm = 7.5;%/(tspan/100); %период повторения модуляций синуса, равный полупериоду меандра
    kamp = Amax/Tm; %коэффициент наклона синуса
    Aspan = mod(kamp*tspan,Amax); %остаток от деления
    
    x(i) = Aspan.*cos(2*pi*tspan); %модулированный синус
    
    Asq = 2; % амплитуда меандра
    y(i) = Asq*square(tspan*pi/Tm); % сам меандр
    
    z(i) = x(i) + y(i); %нужный сигнал/chirp-сигнал
end

% z = transpose(z); %нужный сигнал/chirp-сигнал
CHIR = z/max(abs(z));
figure('Name','Chirp')
plot(n/1000,CHIR)
xlabel('[n] \times 1000, отсчетов')
ylabel('Нормированная амплитуда')
set(gca,'FontName','Times New Roman','FontSize',20)
grid on
%% %%%%%%%%%%%% Генерация ШПС и СШПС с использованием функции chirp %%%%%

if hL == 0.005
    nnn = 100;
else nnn = 3000;
end

LFM2 = chirp(n/nnn, w*2, tendL, 6*w);%, 'logarithmic'); %ШПС
SWL2 = chirp(n/nnn, w*1.2, tendL, 13*w);%, 'logarithmic'); %СШПС

% figure()
% plot(n, LFM2, n, SWL2)
% xlabel('time')
% ylabel('Amplitude')
% grid on

% y = chirp(t,f0,t1,f1) 
% generates samples of a linear swept-frequency cosine signal 
% at the time instances defined in array t. 
% The instantaneous frequency at time 0 is f0, 
% and the instantaneous frequency at time t1 is f1.

%%%%%%%%%%%%%%%%%%%%% Старый код %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% ЛЧМ сигнал %%%%%%%%%%%%%%%%%%%%%
% waitbar(0.56,f,'Calculating LFM...');
% speedL = length(t)/1.5; %коэффициент уменьшения частоты ЛЧМ (во сколько раз - число в знаменателе)
% for i = 1:tendL/hL+1
%     wL(i) = 5*w/(5 + i/speedL);
% end
% wL = transpose(wL);
% Lf = cos(5*pi.*wL.*t);
% LFM = 0.5*Lf./max(abs(Lf));
% LFM = transpose(LFM);
% % figure(5)
% % plot(t,Lf)
% 
% 
% %%%%%%%%%%%%%%%%%% СШЛЧМ сигнал %%%%%%%%%%%%%%%%%%%%%
% % waitbar(0.70,f,'Calculating SWLFM...');
% speedL2 = length(t)/4; %коэффициент уменьшения частоты СШЛЧМ (во сколько раз - число в знаменателе)
% for i = 1:tendL/hL+1
%     wL2(i) = 5*w/(5 + i/speedL2);
% end
% wL2 = transpose(wL2);
% SWL = cos(5*pi.*wL2.*t);
% SWL = 0.5*SWL./max(abs(SWL));
% SWL = transpose(SWL);
% % figure(6)
% plot(t,SWL)

%% %%%%%%%%%%% Гармонический сигнал, пусть косинус %%%%%%%%%%%%%%%%%%%%

coss = cos(n*hL*2*pi);

%%       %%%%%%%%%%%%%%%%%%%%%% PLOTS OF SIGNALS %%%%%%%%%%%%%%%%%%%%
% LFM = LFM2;
% SWL = SWL2;

sh = -2;
signals = [PLorenz + 0 * sh;...
    SprottA   + 1 * sh;...
    SprottB   + 2 * sh;...
    SprottC   + 3 * sh;...
    SprottE   + 4 * sh;...
    Four_wing + 5 * sh;...
    CHIR      + 6 * sh;...
%     %LFM       + 6 * sh;...
    LFM2      + 7 * sh;...
%     %SWL       + 7 * sh;...
    SWL2      + 8 * sh;...
    coss      + 9 * sh;];

figure('Name','Signals')
plot(n/1000, signals);
% legend('Прото-Лоренц', 'Sprott B', 'Новый 3d 4wing','Chirp','ЛЧМ','СШЛЧМ', 'Гармонический'), grid;
legend('Прото-Лоренц','Sprott A', 'Sprott B', 'Sprott C', 'Sprott E', 'Новый 3d 4wing','АМС','ЛЧМ','СШЛЧМ', 'Гармонический'), grid;
set(gca, 'Ylim', [-20, 1],'FontName','Times New Roman','FontSize',20)
xlabel('[n] \times 1000, отсчетов')

%% %%%%%%%%%%%% Фильтрация сигналов средой%%%%%%%%%%%%
% Fsamp = 55000; %типо "сэмпловая" частота (?)
Fsamp = Fsamp*1000;

len = length(t);
[Filtered, pf] = BandTransducer(Fsamp, PLorenz, SprottA, SprottB, SprottC, SprottE, Four_wing, CHIR, LFM2, SWL2, coss);
% waitbar(85,f,'Calculating Spectrum...');

figure('Name','Плотность мощности спектра');
plot(pf, Filtered(:,1),'.-',...% 'k',...
     pf, Filtered(:,2),'-',...% 'k',...
     pf, Filtered(:,3),'--',...
     pf, Filtered(:,4),'--',...
     pf, Filtered(:,5),'--',...
     pf, Filtered(:,6),'-',...% 'k',...
     pf, Filtered(:,7),'-',...% 'k',...
     pf, Filtered(:,8),'--',...% 'k',...
     pf, Filtered(:,9),'k',...% 'k',...
     pf, Filtered(:,10), ':b'); %w, mag2db(ppx2),w, mag2db(ppx3)
 
legend('Прото-Лоренц', 'Sprott A', 'Sprott B', 'Sprott C', 'Sprott E', 'Новый 3d 4wing','АМС','ЛЧМ','СШЛЧМ', 'Гармонический');

%верхнее значение - для больших частот и большего времени соответственно
% set(gca, 'XScale', 'log', 'XLim',[0, 1000],'XTick', [10 30 75 110 150 1000], 'YLim',[-250, -50]), grid;
set(gca, 'XScale', 'log', 'XLim',[0, 1000],'XTick', [1 7 15 30 55 70 100 500], 'YLim',[-250, -50],...
    'FontName','Times New Roman','FontSize',20), grid;

ylabel('Плотность мощности спектра, дБ')
xlabel('Частота, кГц')


%% %%%%%%%%%%%%%%%%%% Моделирование затухания %%%%%%%%%%%%%%%%%%%%%%%%%
% if Fsamp == 30000
%     r = 7000; %в метрах
% else
%     r = 2000; %в метрах
% end
r = r*1000;
T = 4; %температура воды в градусах Цельсия
S = 20; %соленость в промилле (ppt)
Z = 0.05; %глубина в километрах
%%%%%% Диапазон частот для моделирования затухания
FkHz = 10:10:810; %kHz немного подредактированное, чтобы задать края edges
FkHz2 = 10:50:810;
%%%%%%%%%%%%%%%%Выгружаем результаты работы функции
[alpha,aa,signals_f,filt_signals] = absorbation_sound_speed(Fsamp, T, S, Z, r, FkHz, PLorenz, SprottA, SprottB, SprottC, SprottE, Four_wing, CHIR, LFM2, SWL2, coss);

%Потери на поглощение при передаче на расстояние r
TL = 20*log(r) + alpha*r*10^(-3);
Mf = alpha*r*10^(-3); %ослабление в дБ

freq = Fsamp/1000;
alpha2 = alpha(1:5:81);
figure('Name','Alpha')
loglog(FkHz,alpha,FkHz2, alpha2,'k*')
% hold on
% loglog(FkHz2, alpha2,'k*')
set(gca, 'FontName','Times New Roman','FontSize',20)
% title('График зависимости коэффициента поглощения от частоты')
% if r == 2000
%     title('(a)')
% else
%     title('(б)')
% end
xlabel('Частота, кГц')
ylabel('\alpha, дБ/км')
grid on

%     figure(11)
%     plot(FkHz, mag2db(aa), '-*k')
%     set(gca, 'XScale', 'log','FontName','Times New Roman','FontSize',12)
% %     title(['Ослабление в зависимости от частоты на расстоянии ', num2str(r/1000),' км'])
%     if r == 2000
%         title('(a)')
%     else
%         title('(б)')
%     end
%     xlabel('\itf\rm, kHz');
%     ylabel('Ослабление, дБ')
%     set(gca, 'XTick', [1 10 30 55 100 300 400], 'XLim', [5 500]);
%     grid on

figure(411); plot(n/1000, signals_f);
title(['Частота: ',num2str(freq),' кГц;','  Расстояние:',num2str(r/1000),' км'])
legend('Прото-Лоренц', 'Sprott A', 'Sprott B', 'Sprott C', 'Sprott E',...
    'Новый 3d 4wing','АМС','ЛЧМ','СШЛЧМ', 'Гармонический'), grid;
set(gca, 'Ylim', [-20, 1],'FontName','Times New Roman','FontSize',20)
xlabel('[n] \times 1000, отсчетов')

%% %%%%%%%% Далее надо добавить шум и проверить на устойчивость к нему

SA = mean(abs(filt_signals),2);
SA2 = max(abs(filt_signals)')';


        figure('Name','Соотношение среднее значение - максимум');
        c = categorical({'Прото-Лоренц', 'Sprott A', 'Sprott B',...
            'Sprott C','Sprott E', 'Новый 3d 4wing','АМС','ЛЧМ','СШЛЧМ'});%, grid;
        bar(c, [SA, SA2]);
        set(gca,'FontName','Times New Roman','FontSize',20), grid;
        legend('Среднее значение', 'Максимум')
        xlabel('[n] \times 1000, отсчетов')

[Err2,SNR_s,bb1,lags,lagz]= calc_dist_noise(filt_signals);
figure('Name','Погрешность распознавания сигнала от соотношения сигнал-шум');
plot(SNR_s(1,(3:16)), Err2(1,(3:16)),'.-',...% 'k',...
     SNR_s(2,(3:16)), Err2(2,(3:16)),'--',...% 'k',...
     SNR_s(3,(3:16)), Err2(3,(3:16)),'-*',...
     SNR_s(4,(3:16)), Err2(4,(3:16)),'--o',...
     SNR_s(5,(3:16)), Err2(5,(3:16)),'-x',...
     SNR_s(6,(3:16)), Err2(6,(3:16)),'-o',...% 'k',...
     SNR_s(7,(3:16)), Err2(7,(3:16)),'-^' ,...% 'k',...][ 
     SNR_s(8,(3:16)), Err2(8,(3:16)),'--',...% 'k',...
     SNR_s(9,(3:16)), Err2(9,(3:16)),'--^')%,...% 'k',...
%      SNR_s(10,(5:16)), Err2(10,(5:16)), ':b');
 legend('Прото-Лоренц', 'Sprott A', 'Sprott B',...
     'Sprott C', 'Sprott E', 'Новый 3d 4wing',...
     'АМС','ЛЧМ','СШЛЧМ'), grid;

set(gca, 'YScale', 'log', 'XScale', 'log','FontName','Times New Roman','FontSize',20);
xlabel('Соотношение "сигнал/шум"');
ylabel('Погрешность, отсчетов')
% 
% figure(778)
% plot(lags,bb1(1,:),...% 'k',...
%      lags,bb1(2,:),'--',...% 'k',...
%      lags,bb1(3,:),'-.',...
%      lags,bb1(4,:),'-',...
%      lags,bb1(5,:),'-.',...
%      lags,bb1(6,:),'-',...% 'k',...
%      lags,bb1(7,:),'-' ,...% 'k',...][ 
%      lags,bb1(8,:),'--',...% 'k',...
%      lags,bb1(9,:),':')%,...% 'k',...
% %      bb1(10,:), lags, ':b');
%  legend('Прото-Лоренц', 'Sprott A', 'Sprott B',...
%      'Sprott C', 'Sprott E', 'Новый 3d 4wing',...
%      'Chirp','ЛЧМ','СШЛЧМ'), grid;
% % title('Кросс-корреляция')
% set(gca,'FontName','Times New Roman','FontSize',16); % 'YScale', 'log', 'XScale', 'log',
% xlabel('Lag, samples')
% ylabel('Cross-correlation')

% 
% lag0 = 0e3;
% lag1 = 1e3;
% lag2 = 2e3;
% lag3 = 3e3;
% lag4 = 4e3;
% lag5 = 5e3;
% lag6 = 6e3;
% lag7 = 7e3;
% lag8 = 8e3;
% 
% b1 = bb1(1,:);
% b2 = bb1(2,:);
% b3 = bb1(3,:);
% b4 = bb1(4,:);
% b5 = bb1(5,:);
% b6 = bb1(6,:);
% b7 = bb1(7,:);
% b8 = bb1(8,:);
% b9 = bb1(9,:);
% % lags = lags/1000;
% 
% figure(778)
% WINW = 15000;
% WINSHIFT = 500;
% h = plot(lags, b1, lags, b2, lags, b3, lags, b4, lags, b5, lags, b6, lags, b7, lags, b8, lags, b9);
% for i=1:9
%     h(i).Color = [0.8 0.8 0.8];
% end
% hold on
% plot(lags((lagz+lag0-WINW+WINSHIFT):(lagz+lag0+WINW+WINSHIFT)), b1((lagz+lag0-WINW+WINSHIFT):(lagz+lag0+WINW+WINSHIFT)));%,...
% plot(lags((lagz+lag1-WINW+WINSHIFT):(lagz+lag1+WINW+WINSHIFT)), b2((lagz+lag1-WINW+WINSHIFT):(lagz+lag1+WINW+WINSHIFT)));
% plot(lags((lagz+lag2-WINW+WINSHIFT):(lagz+lag2+WINW+WINSHIFT)), b3((lagz+lag2-WINW+WINSHIFT):(lagz+lag2+WINW+WINSHIFT)));
% plot(lags((lagz+lag3-WINW+WINSHIFT):(lagz+lag3+WINW+WINSHIFT)), b4((lagz+lag3-WINW+WINSHIFT):(lagz+lag3+WINW+WINSHIFT)));
% plot(lags((lagz+lag4-WINW+WINSHIFT):(lagz+lag4+WINW+WINSHIFT)), b5((lagz+lag4-WINW+WINSHIFT):(lagz+lag4+WINW+WINSHIFT)));
% plot(lags((lagz+lag5-WINW+WINSHIFT):(lagz+lag5+WINW+WINSHIFT)), b6((lagz+lag5-WINW+WINSHIFT):(lagz+lag5+WINW+WINSHIFT)));
% plot(lags((lagz+lag6-WINW+WINSHIFT):(lagz+lag6+WINW+WINSHIFT)), b7((lagz+lag6-WINW+WINSHIFT):(lagz+lag6+WINW+WINSHIFT)));
% plot(lags((lagz+lag7-WINW+WINSHIFT):(lagz+lag7+WINW+WINSHIFT)), b8((lagz+lag7-WINW+WINSHIFT):(lagz+lag7+WINW+WINSHIFT)));
% plot(lags((lagz+lag8-WINW+WINSHIFT):(lagz+lag8+WINW+WINSHIFT)), b9((lagz+lag8-WINW+WINSHIFT):(lagz+lag8+WINW+WINSHIFT)));
% 
% legend('Прото-Лоренц', 'Sprott A', 'Sprott B',...
%      'Sprott C', 'Sprott E', 'Новый 3d 4wing',...
%      'АМС','ЛЧМ','СШЛЧМ'), grid;
% % title('Кросс-корреляция')
% set(gca,'XLim', [-2000 5000],'FontName','Times New Roman','FontSize',20); % 'YScale', 'log', 'XScale', 'log',
% xlabel('Lag, samples')
% ylabel('Cross-correlation')


















