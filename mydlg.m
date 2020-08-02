function [Fsamp,r,FLAG]=mydlg(xx1,xx2,flag)
% основная функция, в которой создается окно и элементы управления
% создание окна, запись указателя на него в переменную hF
hF=figure('MenuBar','none', 'Position',[150 150 250 200],'NumberTitle','off',...
'Name','Parameters');

% создание статического текста 'Частота='
uicontrol('Style','text','Position',[2 120 140 25],'String','Частота, кГц ',...
    'FontName','Times New Roman', 'FontSize', 14)

% создание области ввода для 'Fsamp', запись в нее значения входного аргумента x1
% и  запись указателя на нее в переменную hInputX1
hInputFsamp=uicontrol('Style','edit','Position',[140 120 60 25],...
'String',num2str(xx1),'FontName','Times New Roman', 'FontSize', 14);

% создание статического текста 'Расстояние ='
uicontrol('Style','text','Position',[5 80 140 25],'String','Расстояние, км',...
    'FontName','Times New Roman', 'FontSize', 14)

% создание области ввода для 'r', запись в нее значения входного аргумента x2
% и  запись указателя на нее в переменную hInputX2
hInputr=uicontrol('Style','edit','Position',[140 80 60 25],...
'String',num2str(xx2),'FontName','Times New Roman', 'FontSize', 14);

% создание флага, запись указателя на него в переменную hFlag
hFlag=uicontrol('Style','checkbox','Position',[5 30 80 15],...
    'String','optimize');
% установка или сброс флага в соответствии со значением входного аргумента flag
set(hFlag,'Value',flag)

% создание кнопки OK, запись указателя на нее в переменную hOK и связь
% ее события Callback со вложенной функцией OK_Callback
hOK=uicontrol('Style','pushbutton','Position',[95 5 70 25],...
'String','OK','FontName','Times New Roman', 'FontSize', 14,'Callback',@OK_Callback);

% создание кнопки Cancel, запись указателя на нее в переменную hCancel и связь
% ее события Callback со вложенной функцией Cancel_Callback
        % hCancel=uicontrol('Style','pushbutton','Position',[175 5 70 25],...
        % 'String','Cancel','Callback',@Cancel_Callback);
        
% приостановка работы основной функции, пока существует диалоговое окно
uiwait(hF)

function OK_Callback(src,evt)
  % вложенная функция для обработки нажатия на кнопку OK

  % считывание значений из областей ввода, преобразование их в числа
  % и запись числовых значений  в выходные аргументы X1 и X2 функции mydlg
  Fsamp = str2num(get(hInputFsamp,'String'));
  r = str2num(get(hInputr,'String'));

    % считывание состояния флага и запись его в выходной аргумент FLAG функции mydlg
    FLAG=get(hFlag,'Value');
    % удаление диалогового окна
    delete(hF)
    end
end
