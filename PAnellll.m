clear all
% hoice = questdlg('А вы точно уверены', 'Warning', 'да', 'нет', 'возможно','');
% f = figure(1);
% panel_1 = uipanel(f,'Title','Main Panel','TitlePosition','centertop','FontSize',12,...
%              'BackgroundColor','white','BorderType','beveledin',...
%              'Position',[.01 .01 0.99 0.99]);
%          
% hbsp = uicontrol('Parent',panel_1,'String','Start calculating',...
%               'Position',[1 1 90 40]);
          
function [X1,X2,FLAG]=mydlg(x1,x2,flag)
% основная функция, в которой создается окно и элементы управления

% создание окна, запись указателя на него в переменную hF
hF=figure('MenuBar','none', 'Position',[150 150 150 100],'NumberTitle','off',...
   	 'Name','Parameters')
              % создание статического текста 'X1=' 
uicontrol('Style','text','Position',[5 70 30 15],'String','X1=')
% создание области ввода для x1, запись в нее значения входного аргумента x1
% и  запись указателя на нее в переменную hInputX1
hInputX1=uicontrol('Style','edit','Position',[40 70 60 15],...
   	 'String',num2str(x1));
              % создание статического текста 'X2='              
uicontrol('Style','text','Position',[5 50 30 15],'String','X2=')
% создание области ввода для x2, запись в нее значения входного аргумента x2
% и  запись указателя на нее в переменную hInputX2
hInputX2=uicontrol('Style','edit','Position',[40 50 60 15],...
    'String',num2str(x2));
% создание флага, запись указателя на него в переменную hFlag
hFlag=uicontrol('Style','checkbox','Position',[5 30 80 15],...
    'String','optimize');
% установка или сброс флага в соответствии со значением входного аргумента flag
set(hFlag,'Value',flag)
% создание кнопки OK, запись указателя на нее в переменную hOK и связь 
% ее события Callback со вложенной функцией OK_Callback
hOK=uicontrol('Style','pushbutton','Position',[5 5 50 15],...
    'String','OK','Callback',@OK_Callback);
% создание кнопки Cancel, запись указателя на нее в переменную hCancel и связь 
% ее события Callback со вложенной функцией Cancel_Callback
hCancel=uicontrol('Style','pushbutton','Position',[65 5 50 15],...
    'String','Cancel','Callback',@Cancel_Callback);
% приостановка работы основной функции, пока существует диалоговое окно
uiwait(hF)
    
    function OK_Callback(src,evt)
    % вложенная функция для обработки нажатия на кнопку OK

    % считывание значений из областей ввода, преобразование их в числа 
    % и запись числовых значений  в выходные аргументы X1 и X2 функции mydlg
    X1=str2num(get(hInputX1,'String'));
    X2=str2num(get(hInputX2,'String'));
    % считывание состояния флага и запись его в выходной аргумент FLAG функции mydlg
    FLAG=get(hFlag,'Value');
    % удаление диалогового окна
    delete(hF)
    end

    function Cancel_Callback(src,evt)
    % вложенная функция для обработки нажатия на кнопку Cancel

   % отказ от всех изменений, записываем в выходные аргументы функции mydlg то, 
   % что было в ее входных аргументах
    X1=x1; X2=x2; FLAG=flag;
    % удаление диалогового окна
    delete(hF)
    end
end



