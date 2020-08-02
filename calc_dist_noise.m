function [Err2,SNR_s,bb1,lags,lagz]= calc_dist_noise(filt_signals)
%%%%%%%%%%%%%% Добавим шум и проверим на устойчивость к нему %%%%%%%%

% lag1 = 0e3;
% lag2 = 2e3;
% lag3 = 4e3;
% lag4 = 3e3;
% lag5 = 1e3;
lagz = 15e3;

SNRs = logspace(-4, 5, 20);

Nexp = 10;
errs = zeros(Nexp, length(SNRs));

L = length(filt_signals(1,:));

koef = L/25001;
lagz = round(koef*15e3);

for i=1:size(filt_signals,1)
    signals_test = filt_signals(i,:);
    for j=1:length(SNRs)
        for k = 1:Nexp
            lag2 = round(rand(1,1)*koef*4e3);
            stub = [zeros(1, lag2), signals_test((1):(L-lag2))] + rand(1,L)./SNRs(j);
            [b1, lags] = crosscorr(signals_test, stub, lagz); [m, id1] = max(b1); %[m, id1] = max(abs(b1));
            errs(k, j) = abs(id1 - lag2 - lagz);
        end
        errs2 = mean(errs, 1);
    end
    Err2(i,:) = errs2; 
    SNR_s(i,:) = SNRs;
    bb1(i,:) = b1;
%     plot(SNRs, errs2, 'LineWidth', 2);
end
end

% figure('Name','Погрешность распознавания сигнала от соотношения сигнал-шум');
% plot(SNR_s(1,(3:16)), Err2(1,(3:16)),'.-',...% 'k',...
%      SNR_s(2,(3:16)), Err2(2,(3:16)),'-',...% 'k',...
%      SNR_s(3,(3:16)), Err2(3,(3:16)),'--',...
%      SNR_s(4,(3:16)), Err2(4,(3:16)),'--',...
%      SNR_s(5,(3:16)), Err2(5,(3:16)),'--',...
%      SNR_s(6,(3:16)), Err2(6,(3:16)),'-',...% 'k',...
%      SNR_s(7,(3:16)), Err2(7,(3:16)),'-',...% 'k',...][ 
%      SNR_s(8,(3:16)), Err2(8,(3:16)),'--',...% 'k',...
%      SNR_s(9,(3:16)), Err2(9,(3:16)),'k')%,...% 'k',...
% %      SNR_s(10,(5:16)), Err2(10,(5:16)), ':b');




%  legend('Прото-Лоренц', 'Sprott A', 'Sprott B',...
%      'Sprott C', 'Sprott E', 'Новый 3d 4wing',...
%      'Chirp','ЛЧМ','СШЛЧМ'), grid;
% 
% set(gca, 'YScale', 'log', 'XScale', 'log');
% xlabel('Соотношение "сигнал/шум"');
% ylabel('Погрешность, отсчетов')
%  
%  
 
 
 
 
 
 
