clear all
% hoice = questdlg('� �� ����� �������', 'Warning', '��', '���', '��������','');
% f = figure(1);
% panel_1 = uipanel(f,'Title','Main Panel','TitlePosition','centertop','FontSize',12,...
%              'BackgroundColor','white','BorderType','beveledin',...
%              'Position',[.01 .01 0.99 0.99]);
%          
% hbsp = uicontrol('Parent',panel_1,'String','Start calculating',...
%               'Position',[1 1 90 40]);
          
function [X1,X2,FLAG]=mydlg(x1,x2,flag)
% �������� �������, � ������� ��������� ���� � �������� ����������

% �������� ����, ������ ��������� �� ���� � ���������� hF
hF=figure('MenuBar','none', 'Position',[150 150 150 100],'NumberTitle','off',...
   	 'Name','Parameters')
              % �������� ������������ ������ 'X1=' 
uicontrol('Style','text','Position',[5 70 30 15],'String','X1=')
% �������� ������� ����� ��� x1, ������ � ��� �������� �������� ��������� x1
% �  ������ ��������� �� ��� � ���������� hInputX1
hInputX1=uicontrol('Style','edit','Position',[40 70 60 15],...
   	 'String',num2str(x1));
              % �������� ������������ ������ 'X2='              
uicontrol('Style','text','Position',[5 50 30 15],'String','X2=')
% �������� ������� ����� ��� x2, ������ � ��� �������� �������� ��������� x2
% �  ������ ��������� �� ��� � ���������� hInputX2
hInputX2=uicontrol('Style','edit','Position',[40 50 60 15],...
    'String',num2str(x2));
% �������� �����, ������ ��������� �� ���� � ���������� hFlag
hFlag=uicontrol('Style','checkbox','Position',[5 30 80 15],...
    'String','optimize');
% ��������� ��� ����� ����� � ������������ �� ��������� �������� ��������� flag
set(hFlag,'Value',flag)
% �������� ������ OK, ������ ��������� �� ��� � ���������� hOK � ����� 
% �� ������� Callback �� ��������� �������� OK_Callback
hOK=uicontrol('Style','pushbutton','Position',[5 5 50 15],...
    'String','OK','Callback',@OK_Callback);
% �������� ������ Cancel, ������ ��������� �� ��� � ���������� hCancel � ����� 
% �� ������� Callback �� ��������� �������� Cancel_Callback
hCancel=uicontrol('Style','pushbutton','Position',[65 5 50 15],...
    'String','Cancel','Callback',@Cancel_Callback);
% ������������ ������ �������� �������, ���� ���������� ���������� ����
uiwait(hF)
    
    function OK_Callback(src,evt)
    % ��������� ������� ��� ��������� ������� �� ������ OK

    % ���������� �������� �� �������� �����, �������������� �� � ����� 
    % � ������ �������� ��������  � �������� ��������� X1 � X2 ������� mydlg
    X1=str2num(get(hInputX1,'String'));
    X2=str2num(get(hInputX2,'String'));
    % ���������� ��������� ����� � ������ ��� � �������� �������� FLAG ������� mydlg
    FLAG=get(hFlag,'Value');
    % �������� ����������� ����
    delete(hF)
    end

    function Cancel_Callback(src,evt)
    % ��������� ������� ��� ��������� ������� �� ������ Cancel

   % ����� �� ���� ���������, ���������� � �������� ��������� ������� mydlg ��, 
   % ��� ���� � �� ������� ����������
    X1=x1; X2=x2; FLAG=flag;
    % �������� ����������� ����
    delete(hF)
    end
end



