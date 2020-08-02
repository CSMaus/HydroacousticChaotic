function [Filtered, pf] = BandTransducer(Fsamp, PLorenz, SprottA, SprottB, SprottC, SprottE, Four_wing, chirp, LFM2, SWL2, coss)

Fs = Fsamp*1000;
lens = Fs;

%%% ��������������� ��������� D ��� ������ tendL � hL

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


% % ������������ ������� ���������� ��������� �������:
% % 
% % Ap - ���������� ��������� � ������ �����������. ����� ���������� Apass.
% % Ast1 - ��������� � ������ ����-��������� � ��������� (������� �� ���������). ����� ���������� Astop1.
% % Ast2 - ��������� �� ������ ����-��������� � ��������� (������� �� ���������). ����� ���������� Astop2.
% % BWp - ������ ������ ����������� �������. ����������� � ������������� �������� �������.
% % BWst - ������ ������ ����������� �������. ����������� � ������������� �������� �������.
% % C - ������������ ���� ������. ��� ��������� ��� ������� ��������� � ������ ����������� ��� ��������� � ������ ������������ ��� ���� �������������� ������� � ����� ��� ���� �� ���� �����.
% % � ������������ �N, Fst1, Fp1, Fp2, Fst2, C� ������ ��������� ����������� � ������� �������� � ������� ����������� ������������. �� ������ ������� ����������� � ����� ��� ���� �������.
% % F3dB1 - ������� ����� ��� ����� �� 3 �� ���� �������� ������ ����������� ��� ������� �����. ����������� � ������������� �������� �������. (���-�������)
% % F3dB2 - ������� ����� ��� ����� �� 3 �� ���� �������� ������ ����������� ��� ������� �����. ����������� � ������������� �������� �������. (���-�������)
% % Fc1 - ������� ����� ��� ����� �� 6 �� ���� �������� ������ ����������� ��� ������� �����. ����������� � ������������� �������� �������. (���-�������)
% % Fc2 - ������� ����� ��� ����� �� 6 �� ���� �������� ������ ����������� ��� ������� �����. ����������� � ������������� �������� �������. (���-�������)
% % Fp1 - ������� �� ���� ������ ������ �����������. ����������� � ������������� �������� �������. ����� ���������� Fpass1.
% % Fp2 - ������� �� ���� ����� ������ �����������. ����������� � ������������� �������� �������. ����� ���������� Fpass2.
% % Fst1 - ������� �� ���� ������ ������� ����-���������. ����������� � ������������� �������� �������. ����� ���������� Fstop1.
% % Fst2 - ������� �� ���� ������ ������� ����-���������. ����������� � ������������� �������� �������. ����� ���������� Fstop2.
% % N - ������� �������� ��� ���-��������. ��� ������� ��������� � ����������� ��� ���-��������, ����� �������� na � nb �� �������.
% % Na - ����������� ������� ��� ���-��������
% % Nb - ������� ��������� ��� ���-��������

pf = pfCh/1000;
Filtered = [ppLor, ppSprA, ppSprB, ppSprC, ppSprE pp4wing, ppChirp, ppLFM, ppSWL, ppcoss];
end