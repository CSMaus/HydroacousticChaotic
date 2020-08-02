function [Fsamp,r,FLAG]=mydlg(xx1,xx2,flag)
% �������� �������, � ������� ��������� ���� � �������� ����������
% �������� ����, ������ ��������� �� ���� � ���������� hF
hF=figure('MenuBar','none', 'Position',[150 150 250 200],'NumberTitle','off',...
'Name','Parameters');

% �������� ������������ ������ '�������='
uicontrol('Style','text','Position',[2 120 140 25],'String','�������, ��� ',...
    'FontName','Times New Roman', 'FontSize', 14)

% �������� ������� ����� ��� 'Fsamp', ������ � ��� �������� �������� ��������� x1
% �  ������ ��������� �� ��� � ���������� hInputX1
hInputFsamp=uicontrol('Style','edit','Position',[140 120 60 25],...
'String',num2str(xx1),'FontName','Times New Roman', 'FontSize', 14);

% �������� ������������ ������ '���������� ='
uicontrol('Style','text','Position',[5 80 140 25],'String','����������, ��',...
    'FontName','Times New Roman', 'FontSize', 14)

% �������� ������� ����� ��� 'r', ������ � ��� �������� �������� ��������� x2
% �  ������ ��������� �� ��� � ���������� hInputX2
hInputr=uicontrol('Style','edit','Position',[140 80 60 25],...
'String',num2str(xx2),'FontName','Times New Roman', 'FontSize', 14);

% �������� �����, ������ ��������� �� ���� � ���������� hFlag
hFlag=uicontrol('Style','checkbox','Position',[5 30 80 15],...
    'String','optimize');
% ��������� ��� ����� ����� � ������������ �� ��������� �������� ��������� flag
set(hFlag,'Value',flag)

% �������� ������ OK, ������ ��������� �� ��� � ���������� hOK � �����
% �� ������� Callback �� ��������� �������� OK_Callback
hOK=uicontrol('Style','pushbutton','Position',[95 5 70 25],...
'String','OK','FontName','Times New Roman', 'FontSize', 14,'Callback',@OK_Callback);

% �������� ������ Cancel, ������ ��������� �� ��� � ���������� hCancel � �����
% �� ������� Callback �� ��������� �������� Cancel_Callback
        % hCancel=uicontrol('Style','pushbutton','Position',[175 5 70 25],...
        % 'String','Cancel','Callback',@Cancel_Callback);
        
% ������������ ������ �������� �������, ���� ���������� ���������� ����
uiwait(hF)

function OK_Callback(src,evt)
  % ��������� ������� ��� ��������� ������� �� ������ OK

  % ���������� �������� �� �������� �����, �������������� �� � �����
  % � ������ �������� ��������  � �������� ��������� X1 � X2 ������� mydlg
  Fsamp = str2num(get(hInputFsamp,'String'));
  r = str2num(get(hInputr,'String'));

    % ���������� ��������� ����� � ������ ��� � �������� �������� FLAG ������� mydlg
    FLAG=get(hFlag,'Value');
    % �������� ����������� ����
    delete(hF)
    end
end
