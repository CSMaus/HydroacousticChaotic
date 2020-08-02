function [Filtered, pf] = BandTransducer(Fsamp, PLorenz, SprottA, SprottB, SprottC, SprottE, Four_wing, chirp, LFM2, SWL2, coss)

Fs = Fsamp*1000;
lens = Fs;

%%% Соответствующее выражение D для других tendL и hL

D = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
    50/Fsamp,55/Fsamp,125/Fsamp,130/Fsamp,2,1,2);
D2 = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
    150/Fsamp,20/Fsamp,70/Fsamp,80/Fsamp,2,1,2);


% D = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
%     10/Fsamp,15/Fsamp,75/Fsamp,80/Fsamp,2,1,2);


BandTrans = design(D,'butter');

F_PLorenz = filter(BandTrans,PLorenz);
% F_PLorenz = PLorenz;
[ppLor,pfL] = pwelch(F_PLorenz,[],[],[],lens);
ppLor = mag2db(ppLor);

F_SprottA = filter(BandTrans,SprottA);
% F_SprottA = SprottA;
[ppSprA,pfS] = pwelch(F_SprottA,[],[],[],lens);
ppSprA = mag2db(ppSprA);

F_SprottB = filter(BandTrans,SprottB);
% F_SprottB = SprottB;
[ppSprB,pfS] = pwelch(F_SprottB,[],[],[],lens);
ppSprB = mag2db(ppSprB);

F_SprottC = filter(BandTrans,SprottC);
% F_SprottC = SprottC;
[ppSprC,pfS] = pwelch(F_SprottC,[],[],[],lens);
ppSprC = mag2db(ppSprC);

F_SprottE = filter(BandTrans,SprottE);
% F_SprottE = SprottE;
[ppSprE,pfS] = pwelch(F_SprottE,[],[],[],lens);
ppSprE = mag2db(ppSprE);

F_Four_wing = filter(BandTrans,Four_wing);
% F_Four_wing = Four_wing;
[pp4wing,pf4] = pwelch(F_Four_wing,[],[],[],lens);
pp4wing = mag2db(pp4wing);

F_chirp = filter(BandTrans,chirp);
% F_chirp = chirp;
[ppChirp,pfCh] = pwelch(F_chirp,[],[],[],lens);
ppChirp = mag2db(ppChirp);

F_LFM = filter(BandTrans,LFM2);
% F_LFM = LFM2;
[ppLFM,pfLFM] = pwelch(F_LFM,[],[],[],lens);
ppLFM = mag2db(ppLFM);

F_SWL = filter(BandTrans,SWL2);
% F_SWL = SWL2;
[ppSWL,pfSWL] = pwelch(F_SWL,[],[],[],lens);
ppSWL = mag2db(ppSWL);

F_coss = filter(BandTrans,coss);
% F_coss = coss;
[ppcoss,pfcoss] = pwelch(F_coss,[],[],[],lens);
ppcoss = mag2db(ppcoss);


% % Спецификации фильтра определены следующим образом:
% % 
% % Ap - допустимая пульсация в полосе пропускания. Также называется Apass.
% % Ast1 - затухание в первом стоп-диапазоне в децибелах (единицы по умолчанию). Также называется Astop1.
% % Ast2 - затухание во втором стоп-диапазоне в децибелах (единицы по умолчанию). Также называется Astop2.
% % BWp - ширина полосы пропускания фильтра. Указывается в нормированных единицах частоты.
% % BWst - ширина полосы пропускания фильтра. Указывается в нормированных единицах частоты.
% % C - Ограниченный флаг группы. Это позволяет вам указать пульсацию в полосе пропускания или затухание в полосе задерживания для схем фиксированного порядка в одной или двух из трех полос.
% % В спецификации «N, Fst1, Fp1, Fp2, Fst2, C» нельзя указывать ограничения в полосах задержек и полосах пропускания одновременно. Вы можете указать ограничения в одной или двух полосах.
% % F3dB1 - частота среза для точки на 3 дБ ниже значения полосы пропускания для первого среза. Указывается в нормированных единицах частоты. (БИХ-фильтры)
% % F3dB2 - частота среза для точки на 3 дБ ниже значения полосы пропускания для второго среза. Указывается в нормированных единицах частоты. (БИХ-фильтры)
% % Fc1 - частота среза для точки на 6 дБ ниже значения полосы пропускания для первого среза. Указывается в нормированных единицах частоты. (КИХ-фильтры)
% % Fc2 - частота среза для точки на 6 дБ ниже значения полосы пропускания для второго среза. Указывается в нормированных единицах частоты. (КИХ-фильтры)
% % Fp1 - частота на краю начала полосы пропускания. Указывается в нормированных единицах частоты. Также называется Fpass1.
% % Fp2 - частота на краю конца полосы пропускания. Указывается в нормированных единицах частоты. Также называется Fpass2.
% % Fst1 - частота на краю начала первого стоп-диапазона. Указывается в нормированных единицах частоты. Также называется Fstop1.
% % Fst2 - частота на краю начала второго стоп-диапазона. Указывается в нормированных единицах частоты. Также называется Fstop2.
% % N - порядок фильтров для КИХ-фильтров. Или порядки числителя и знаменателя для БИХ-фильтров, когда значения na и nb не указаны.
% % Na - знаменатель порядка для БИХ-фильтров
% % Nb - порядок нумерации для БИХ-фильтров

pf = pfCh/1000;
Filtered = [ppLor, ppSprA, ppSprB, ppSprC, ppSprE pp4wing, ppChirp, ppLFM, ppSWL, ppcoss];
end