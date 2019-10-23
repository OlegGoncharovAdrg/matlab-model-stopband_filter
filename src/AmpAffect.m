

function varargout = AmpAffect(varargin)
% AMPAFFECT M-file for AmpAffect.fig
%      AMPAFFECT, by itself, creates a new AMPAFFECT or raises the existing
%      singleton*.
%
%      H = AMPAFFECT returns the handle to a new AMPAFFECT or the handle to
%      the existing singleton*.
%
%      AMPAFFECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AMPAFFECT.M with the given input arguments.
%
%      AMPAFFECT('Property','Value',...) creates a new AMPAFFECT or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AmpAffect_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AmpAffect_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AmpAffect

% Last Modified by GUIDE v2.5 01-Sep-2010 11:01:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AmpAffect_OpeningFcn, ...
                   'gui_OutputFcn',  @AmpAffect_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before AmpAffect is made visible.
function AmpAffect_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AmpAffect (see VARARGIN)

% Choose default command line output for AmpAffect
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);

global c
c = 3e8;  % �������� �����

% ������� �������
set(handles.panLoadSrcData, 'position', [0.02 0.02 0.96 0.85]);
set(handles.panAcf,         'position', [0.02 0.02 0.96 0.85]);
set(handles.panRcvParam,    'position', [0.02 0.02 0.96 0.85]);
set(handles.panRcvMdl,      'position', [0.02 0.02 0.96 0.85]);
set(handles.panRfAffect,    'position', [0.02 0.02 0.96 0.85]);

% ���������� ����������
set(handles.panAcfProgr,     'visible', 'off');
set(handles.lblAcfProgr,     'position', [0 0.5 0.01 0.45]);
set(handles.lblAcfProgrDscr, 'position', [0 0 1 0.45]);

set(handles.panRcvMdlProgr,     'visible', 'off');
set(handles.lblRcvMdlProgr,     'position', [0 0.5 0.01 0.45]);
set(handles.lblRcvMdlProgrDscr, 'position', [0 0 1 0.45]);


% ������ ������ - �������� ����� ������. ������ ������ �������� �� ��������� ��� �������.
set(handles.panLoadSrcData, 'visible', 'on');

set(handles.tabSrcData, 'value', 1);

set(handles.tabAcf, 'value', 0);
set(handles.panAcf, 'visible', 'off');

set(handles.tabRcvParam, 'value', 0);
set(handles.panRcvParam, 'visible', 'off');

set(handles.tabRcvMdl, 'value', 0);
set(handles.panRcvMdl, 'visible', 'off');

set(handles.tabRfAffect, 'value', 0);
set(handles.panRfAffect, 'visible', 'off');

% ���� �� ��������� �������� ������, ��������� ������� ���������
set(handles.tabAcf, 'enable', 'off');
set(handles.tabRcvParam, 'enable', 'off');
set(handles.tabRcvMdl, 'enable', 'off');
set(handles.tabRfAffect, 'enable', 'off');


% % �� ��������� ������ �������� L1
% set(handles.buttL1, 'value', 1);
% set(handles.buttL2, 'value', 0);

% �� ��������� ������ ������ ��
set(handles.buttPT, 'value', 0);
set(handles.buttVT, 'value', 1);

set(handles.buttPT, 'string', 'BOC(1, 1)');
set(handles.buttVT, 'string', 'BOC(5, 2.5)');

set(handles.buttPTRf, 'value', 0);
set(handles.buttVTRf, 'value', 1);

set(handles.buttPTRf, 'string', 'BOC(1, 1)');
set(handles.buttVTRf, 'string', 'BOC(5, 2.5)');

set(handles.edRfAffectTime, 'string', '1');


% uiwait makes AmpAffect wait for user response (see UIRESUME)
% uiwait(handles.figure1);


global f f0 K AMP PH
f0 = 1602e6;
f = [];
K = [];
AMP = [];
PH = [];


global UsrRcvParams PredRcv
UsrRcvParams = [10 10 0.5;   % ���� �-822
                10 50 3;   % ��������
                5 30 2;    % ��������� ��������, ����������� �������������
                7 25 1];   % ������ ��� (���), ������ ��� � ��� (��)

PredRcv = 2;               % ������ ������� ��������� 2 ���������, �� ������ ������
RcvList = {'���� �-822' '��������'}; % ��� �� ��������
for i=3:size(UsrRcvParams, 1)
    RcvList{i} = sprintf('���������������� %2d\n', i-PredRcv);
end

set(handles.popRcvType, 'string', RcvList);

% �� ��������� ���������� ��������������� �������, ��������� �������� ������ ������
set(handles.edRcvIntBand, 'enable', 'off');
set(handles.edRcvPllBand, 'enable', 'off');
set(handles.edRcvDllBand, 'enable', 'off');

% � ���� ������������� ����� ���� � ������������ ������ ������
set(handles.edRcvMdlMFreq, 'enable', 'off');
set(handles.edRcvMdlDFreq, 'enable', 'off');
set(handles.edRcvMdlMTau, 'enable', 'off');
set(handles.edRcvMdlDTau, 'enable', 'off');

set(handles.chkIFAoff, 'value', 0);
set(handles.chkIFAoff, 'string', '����');

% ������ ���������� ��������������
set(handles.buttRcvMdlStart, 'Enable', 'on');
set(handles.buttRcvMdlStop, 'Enable', 'off');

% �� ��������� ����� ������������ ���������
set(handles.buttTauMode, 'value', 1);

% ������ �����
set(handles.text62, 'string', 'Frequency');
set(handles.text76, 'string', 'Frequency');
for i=-10:10
    Lits{i+11}=(1565+i)*1.023;
end

set(handles.popLit, 'string', Lits);
set(handles.popLit, 'value', 11);  % �� ��������� ������� ������
popLit_Callback(handles.popLit, [], handles);

set(handles.popLitRf, 'string', Lits);
set(handles.popLitRf, 'value', 11);  % �� ��������� ������� ������
popLit_Callback(handles.popLitRf, [], handles);


global RcvMdlDisables % ������ ���������, ������� ����� ������������� ��� ������ �������������
RcvMdlDisables = [handles.popRcvType;
                  handles.edRcvIntBand;
                  handles.edRcvPllBand;
                  handles.edRcvDllBand;
                  handles.popLit;
                  handles.popLitRf;                  
                 ];

set(handles.buttSaveAcf, 'enable', 'off');
set(handles.buttSaveErrFreq, 'enable', 'off');
set(handles.buttSaveErrTau, 'enable', 'off');


% ������ �����������
set(handles.panRcvParam, 'visible', 'off');
set(handles.panRcvMdl, 'visible', 'off');
set(handles.tabRcvParam, 'visible', 'off');
set(handles.tabRcvMdl, 'visible', 'off');

set(handles.edAcfTauMax, 'String', 300);

SrcDataInit(handles, 'data/afr.csv');
SrcDataInitPFR(handles, 'data/pfr.csv');


% --- Outputs from this function are returned to the command line.
function varargout = AmpAffect_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function density_Callback(hObject, eventdata, handles)
% hObject    handle to density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of density as text
%        str2double(get(hObject,'String')) returns contents of density as a double
density = str2double(get(hObject, 'String'));
if isnan(density)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
%handles.metricdata.density = density;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function volume_CreateFcn(hObject, eventdata, handles)
% hObject    handle to volume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function volume_Callback(hObject, eventdata, handles)
% hObject    handle to volume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of volume as text
%        str2double(get(hObject,'String')) returns contents of volume as a double
volume = str2double(get(hObject, 'String'));
if isnan(volume)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new volume value
handles.metricdata.volume = volume;
guidata(hObject,handles)

% --- Executes on button press in calculate.
function calculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%mass = handles.metricdata.density * handles.metricdata.volume;
set(handles.mass, 'String', mass);

% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

initialize_gui(gcbf, handles, true);

% --- Executes when selected object changed in unitgroup.
function unitgroup_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (hObject == handles.english)
    set(handles.text4, 'String', 'lb/cu.in');
    set(handles.text5, 'String', 'cu.in');
    set(handles.text6, 'String', 'lb');
else
    set(handles.text4, 'String', 'kg/cu.m');
    set(handles.text5, 'String', 'cu.m');
    set(handles.text6, 'String', 'kg');
end

% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.
if isfield(handles, 'metricdata') && ~isreset
    return;
end

%handles.metricdata.density = 0;
handles.metricdata.volume  = 0;

%set(handles.density, 'String', handles.metricdata.density);
%set(handles.volume,  'String', handles.metricdata.volume);
%set(handles.mass, 'String', 0);

%set(handles.unitgroup, 'SelectedObject', handles.english);

%set(handles.text4, 'String', 'lb/cu.in');
%set(handles.text5, 'String', 'cu.in');
%set(handles.text6, 'String', 'lb');

% Update handles structure
guidata(handles.figure1, handles);



% % --- Executes on button press in buttL1.
% function buttL1_Callback(hObject, eventdata, handles)
% % hObject    handle to buttL1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)

% % Hint: get(hObject,'Value') returns toggle state of buttL1
% set(handles.buttL1, 'Value', 1);
% set(handles.buttL2, 'Value', 0);

% global f0 f
% f0 = 1602e6;


% if (size(f, 2)>1)
%     if (f0<min(f)) | (f0>max(f))
%         msgbox(sprintf('��������������: ��������� ��������� �������������� �� ��������� � ��������� L1 (f0=%7.2f ���)', f0/1e6));
%     end
% end

% popLit_Callback(handles.popLit, [], handles);

% --- Executes on button press in tabSrcData.
function tabSrcData_Callback(hObject, eventdata, handles)
% hObject    handle to tabSrcData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    set(handles.tabSrcData, 'value', 1);
    set(handles.tabAcf, 'value', 0);
    set(handles.tabRcvParam, 'value', 0);
    set(handles.tabRcvMdl, 'value', 0);
    set(handles.tabRfAffect, 'value', 0);
    
    set(handles.panLoadSrcData, 'visible', 'on');
    set(handles.panAcf, 'visible', 'off');
    set(handles.panRcvParam, 'visible', 'off');
    set(handles.panRcvMdl, 'visible', 'off');
    set(handles.panRfAffect, 'visible', 'off');        

% --- Executes on button press in tabAcf.
function tabAcf_Callback(hObject, eventdata, handles)
% hObject    handle to tabAcf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    set(handles.tabSrcData, 'value', 0);
    set(handles.tabAcf, 'value', 1);
    set(handles.tabRcvParam, 'value', 0);
    set(handles.tabRcvMdl, 'value', 0);
    set(handles.tabRfAffect, 'value', 0);
    
    set(handles.panLoadSrcData, 'visible', 'off');    
    set(handles.panAcf, 'visible', 'on');
    set(handles.panRcvParam, 'visible', 'off');
    set(handles.panRcvMdl, 'visible', 'off');
    set(handles.panRfAffect, 'visible', 'off');        


% --- Executes on button press in tabRcvParam.
function tabRcvParam_Callback(hObject, eventdata, handles)
% hObject    handle to tabRcvParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    set(handles.tabSrcData, 'value', 0);
    set(handles.tabAcf, 'value', 0);
    set(handles.tabRcvParam, 'value', 1);
    set(handles.tabRcvMdl, 'value', 0);
    set(handles.tabRfAffect, 'value', 0);

    set(handles.panLoadSrcData, 'visible', 'off');
    set(handles.panAcf, 'visible', 'off');
    set(handles.panRcvParam, 'visible', 'on');
    set(handles.panRcvMdl, 'visible', 'off');
    set(handles.panRfAffect, 'visible', 'off');        

    
% --- Executes on button press in buttL2.
% function buttL2_Callback(hObject, eventdata, handles)
% % hObject    handle to buttL2 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)

% % Hint: get(hObject,'Value') returns toggle state of buttL2
% set(handles.buttL1, 'Value', 0);
% set(handles.buttL2, 'Value', 1);

% global f0 f
% f0 = 1247e6;



if (size(f, 2)>1)
    if (f0<min(f)) | (f0>max(f))
        msgbox(sprintf('��������������: ��������� ��������� �������������� �� ��������� � ��������� L2 (f0=%7.2f ���)', f0/1e6));
    end
end

popLit_Callback(handles.popLit, [], handles);
popLitRf_Callback(handles.popLitRf, [], handles);

% --- Executes on button press in buttPT.
function buttPT_Callback(hObject, eventdata, handles)
% hObject    handle to buttPT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of buttPT
    set(handles.buttPT, 'Value', 1);
    set(handles.buttVT, 'Value', 0);
    set(handles.buttPTRf, 'Value', 1);
    set(handles.buttVTRf, 'Value', 0);
    
    set(handles.edAcfTauMax, 'String', 3000);
    

% --- Executes on button press in buttVT.
function buttVT_Callback(hObject, eventdata, handles)
% hObject    handle to buttVT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of buttVT
    set(handles.buttPT, 'Value', 0);
    set(handles.buttVT, 'Value', 1);
    set(handles.buttPTRf, 'Value', 0);
    set(handles.buttVTRf, 'Value', 1);

    set(handles.edAcfTauMax, 'String', 300);

function SrcDataInit(handles, sFileName)

    global f f0 AMP PH c K
    
    fid = fopen(sFileName, 'r');

%    fscanf(fid, '%s\n%s\n%s\n');  % ��� �������� �������
    z = fscanf(fid, '%g, %g, %g\n', [3 1]);
    for i=1:1000
        if size(z, 1) == 3
            break;
        else
            fscanf(fid, '%c', 1);
            z = fscanf(fid, '%g, %g, %g\n', [3 1]);
        end
    end
    
    y = [z fscanf(fid, '%g, %g, %g\n', [3 Inf])];
    fclose(fid);
    f = y(1, :) + 305e3;  % �������� �� ���� ��������� epsilon � �������
    df = diff(f); df = [df(1) df];
    AMP = 10.^(y(2, :)/20);;
    %    K = y(2, :) .* exp(-1i*cumsum( 2*pi*gd .* df ));
    %    K = y(2, :) + 1i*y(3, :);

   
    set(handles.lblSrcData, 'string', sFileName);
    
    % ��������� ������ ��� �������� ������������ �� ������
    plot(handles.pltAFR, f/1e6, 20*log10(abs(AMP)));
    grid(handles.pltAFR, 'on');

    if (size(PH, 2) > 0) & size(AMP, 2) == size(PH, 2)    
        K = AMP .* exp(1i * PH);
        ph = phase(K);
        ax = plot(handles.pltPFR, f/1e6, ph/pi*180);        
    else
        ax = plot(handles.pltPFR, 0, 0);
    end
    grid(handles.pltPFR, 'on');    
    
    
    set(handles.pltAFR, 'ButtonDownFcn',  @pltAFR_ButtonDownFcn);
    
    % ax = plotyy(handles.pltPFR, f/1e6, gd/1e-9, f/1e6, gd*c);
    % ylims = ylim(ax(1));
    % ylim(ax(2), ylims*1e-9*c);
    
    % grid(handles.pltPFR, 'on');

    %    set(handles.pltPFR, 'ButtonDownFcn',  @pltPFR_ButtonDownFcn);
    % ���� ����� ���������� ������� ������ �������, �������, ��� ������ ��������� ���������
    if (size(AMP, 2) > 1) & (size(PH, 2) > 1) & (size(AMP, 2) == size(PH, 2))
        % ����� ���������� ������ ����� ��������� ��������� ������� ���������
        set(handles.tabAcf, 'enable', 'on');
        set(handles.tabRcvParam, 'enable', 'on');
        set(handles.tabRcvMdl, 'enable', 'on');     
        set(handles.tabRfAffect, 'enable', 'on');
        K = AMP .* exp(1i * PH);
    else
        % ���� ������ ������ ��������� ��������, ����� ��������� ������� ���������
        set(handles.tabAcf, 'enable', 'off');
        set(handles.tabRcvParam, 'enable', 'off');
        set(handles.tabRcvMdl, 'enable', 'off'); 
        set(handles.tabRfAffect, 'enable', 'off');        
        K = [];
    end
    
    
% --- Executes on button press in buttLoadSrcData.
function buttLoadSrcData_Callback(hObject, eventdata, handles)
% hObject    handle to buttLoadSrcData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ���������� ����� ���
    
    global f f0 AMP PH c K
    
    % ����������� ��� ����� � ����������� �������
    [sFileName sPathName] = uigetfile('*.csv');

    SrcDataInit(handles, [sPathName sFileName]);
    
    % ���� � ������� CSV, ��� ������� - �������, ���������, ����
    % fid = fopen([sPathName sFileName], 'r');
    
    % %    fscanf(fid, '%s\n%s\n%s\n');  % ��� �������� �������
    % z = fscanf(fid, '%g, %g, %g\n', [3 1]);
    % for i=1:1000
    %     if size(z, 1) == 3
    %         break;
    %     else
    %         fscanf(fid, '%c', 1);
    %         z = fscanf(fid, '%g, %g, %g\n', [3 1]);
    %     end
    % end
    
    % y = [z fscanf(fid, '%g, %g, %g\n', [3 Inf])];
    % fclose(fid);
    % f = y(1, :) + 305e3;  % �������� �� ���� ��������� epsilon � �������
    % df = diff(f); df = [df(1) df];
    % AMP = 10.^(y(2, :)/20);;
    % %    K = y(2, :) .* exp(-1i*cumsum( 2*pi*gd .* df ));
    % %    K = y(2, :) + 1i*y(3, :);

   
    % set(handles.lblSrcData, 'string', sFileName);
    
    % % ��������� ������ ��� �������� ������������ �� ������
    % plot(handles.pltAFR, f/1e6, 20*log10(abs(AMP)));
    % grid(handles.pltAFR, 'on');

    % if (size(PH, 2) > 0) & size(AMP, 2) == size(PH, 2)    
    %     K = AMP .* exp(1i * PH);
    %     ph = phase(K);
    %     ax = plot(handles.pltPFR, f/1e6, ph/pi*180);        
    % else
    %     ax = plot(handles.pltPFR, 0, 0);
    % end
    % grid(handles.pltPFR, 'on');    
    
    
    % set(handles.pltAFR, 'ButtonDownFcn',  @pltAFR_ButtonDownFcn);
    
    % % ax = plotyy(handles.pltPFR, f/1e6, gd/1e-9, f/1e6, gd*c);
    % % ylims = ylim(ax(1));
    % % ylim(ax(2), ylims*1e-9*c);
    
    % % grid(handles.pltPFR, 'on');

    % %    set(handles.pltPFR, 'ButtonDownFcn',  @pltPFR_ButtonDownFcn);
    % % ���� ����� ���������� ������� ������ �������, �������, ��� ������ ��������� ���������
    % if (size(AMP, 2) > 1) & (size(PH, 2) > 1) & (size(AMP, 2) == size(PH, 2))
    %     % ����� ���������� ������ ����� ��������� ��������� ������� ���������
    %     set(handles.tabAcf, 'enable', 'on');
    %     set(handles.tabRcvParam, 'enable', 'on');
    %     set(handles.tabRcvMdl, 'enable', 'on');     
    %     set(handles.tabRfAffect, 'enable', 'on');
    %     K = AMP .* exp(1i * PH);
    % else
    %     % ���� ������ ������ ��������� ��������, ����� ��������� ������� ���������
    %     set(handles.tabAcf, 'enable', 'off');
    %     set(handles.tabRcvParam, 'enable', 'off');
    %     set(handles.tabRcvMdl, 'enable', 'off'); 
    %     set(handles.tabRfAffect, 'enable', 'off');        
    %     K = [];
    % end
    
    % if (size(f, 2)>1)
    %     if (f0<min(f)) | (f0>max(f))  % ��� �� ������� ��������� ��������
    %         if get(handles.buttL1, 'value') == 1  % ������� �������� - L1
    %             buttL2_Callback(handles.buttL1, [], handles); % ������������� � L2
    %             %                msgbox(sprintf('��������������: ��������� ��������� �������������� �� ��������� � ��������� L1 (f0=%7.2f ���)', f0/1e6));
                
    %         else                                  % ������� �������� - L2
    %             buttL1_Callback(handles.buttL1, [], handles); % ������������� � L1
    %             %                msgbox(sprintf('��������������: ��������� ��������� �������������� �� ��������� � ��������� L2 (f0=%7.2f ���)', f0/1e6));
    %         end
    %     end
    % end


function SrcDataInitPFR(handles, sFileName)

    global f f0 AMP PH c K
    
    % ���� � ������� CSV, ��� ������� - �������, ���������, ����
    fid = fopen(sFileName, 'r');
    
    %    fscanf(fid, '%s\n%s\n%s\n');  % ��� �������� �������
    z = fscanf(fid, '%g, %g, %g\n', [3 1]);
    for i=1:1000
        if size(z, 1) == 3
            break;
        else
            fscanf(fid, '%c', 1);
            z = fscanf(fid, '%g, %g, %g\n', [3 1]);
        end
    end
    
    y = [z fscanf(fid, '%g, %g, %g\n', [3 Inf])];
    fclose(fid);
    f = y(1, :) + 305e3;  % �������� �� ���� ��������� epsilon � �������    
    df = diff(f); df = [df(1) df];
    PH = y(2, :)/180*pi;

    %    K = y(2, :) .* exp(-1i*cumsum( 2*pi*gd .* df ));
    %    K = y(2, :) + 1i*y(3, :);

    set(handles.lblSrcDataPfr, 'string', sFileName);
    
    % ��������� ������ ��� �������� ������������ �� ������
    % plot(handles.pltAFR, f/1e6, 20*log10(abs(AMP)));
    % grid(handles.pltAFR, 'on');
    
    % set(handles.pltAFR, 'ButtonDownFcn',  @pltAFR_ButtonDownFcn);
    if (size(AMP, 2) > 0) &  size(AMP, 2) == size(PH, 2)    
        K = AMP .* exp(1i * PH);
        ph = phase(K);
        ax = plot(handles.pltPFR, f/1e6, ph/pi*180);        
    else
        ax = plot(handles.pltPFR, f/1e6, PH/pi*180);
    end
            
    %    ylims = ylim(ax(1));
    %    ylim(ax(2), ylims*1e-9*c);
    
    grid(handles.pltPFR, 'on');

    set(handles.pltPFR, 'ButtonDownFcn',  @pltPFR_ButtonDownFcn);
    
    % ���� ����� ���������� ������� ������ �������, �������, ��� ������ ��������� ���������
    if (size(AMP, 2) > 1) & (size(PH, 2) > 1) & (size(AMP, 2) == size(PH, 2))
        % ����� ���������� ������ ����� ��������� ��������� ������� ���������
        set(handles.tabAcf, 'enable', 'on');
        set(handles.tabRcvParam, 'enable', 'on');
        set(handles.tabRcvMdl, 'enable', 'on');     
        set(handles.tabRfAffect, 'enable', 'on');        
        K = AMP .* exp(1i * PH);
    else
        % ���� ������ ������ ��������� ��������, ����� ��������� ������� ���������
        set(handles.tabAcf, 'enable', 'off');
        set(handles.tabRcvParam, 'enable', 'off');
        set(handles.tabRcvMdl, 'enable', 'off'); 
        set(handles.tabRfAffect, 'enable', 'off');        
        K = [];
    end

    
% --- Executes on button press in buttLoadSrcDataPFR.
function buttLoadSrcDataPFR_Callback(hObject, eventdata, handles)
% hObject    handle to buttLoadSrcDataPFR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ���������� ����� ���
    
    % global f f0 AMP PH c K
    
    % ����������� ��� ����� � ����������� �������
    [sFileName sPathName] = uigetfile('*.csv');

    SrcDataInitPFR(handles, [sFileName sPathName]);

    
    
    % % ���� � ������� CSV, ��� ������� - �������, ���������, ����
    % fid = fopen([sPathName sFileName], 'r');
    
    % %    fscanf(fid, '%s\n%s\n%s\n');  % ��� �������� �������
    % z = fscanf(fid, '%g, %g, %g\n', [3 1]);
    % for i=1:1000
    %     if size(z, 1) == 3
    %         break;
    %     else
    %         fscanf(fid, '%c', 1);
    %         z = fscanf(fid, '%g, %g, %g\n', [3 1]);
    %     end
    % end
    
    % y = [z fscanf(fid, '%g, %g, %g\n', [3 Inf])];
    % fclose(fid);
    % f = y(1, :) + 305e3;  % �������� �� ���� ��������� epsilon � �������    
    % df = diff(f); df = [df(1) df];
    % PH = y(2, :)/180*pi;

    % %    K = y(2, :) .* exp(-1i*cumsum( 2*pi*gd .* df ));
    % %    K = y(2, :) + 1i*y(3, :);

    % set(handles.lblSrcDataPfr, 'string', sFileName);
    
    % % ��������� ������ ��� �������� ������������ �� ������
    % % plot(handles.pltAFR, f/1e6, 20*log10(abs(AMP)));
    % % grid(handles.pltAFR, 'on');
    
    % % set(handles.pltAFR, 'ButtonDownFcn',  @pltAFR_ButtonDownFcn);
    % if (size(AMP, 2) > 0) &  size(AMP, 2) == size(PH, 2)    
    %     K = AMP .* exp(1i * PH);
    %     ph = phase(K);
    %     ax = plot(handles.pltPFR, f/1e6, ph/pi*180);        
    % else
    %     ax = plot(handles.pltPFR, f/1e6, PH/pi*180);
    % end
            
    % %    ylims = ylim(ax(1));
    % %    ylim(ax(2), ylims*1e-9*c);
    
    % grid(handles.pltPFR, 'on');

    % set(handles.pltPFR, 'ButtonDownFcn',  @pltPFR_ButtonDownFcn);
    
    % % ���� ����� ���������� ������� ������ �������, �������, ��� ������ ��������� ���������
    % if (size(AMP, 2) > 1) & (size(PH, 2) > 1) & (size(AMP, 2) == size(PH, 2))
    %     % ����� ���������� ������ ����� ��������� ��������� ������� ���������
    %     set(handles.tabAcf, 'enable', 'on');
    %     set(handles.tabRcvParam, 'enable', 'on');
    %     set(handles.tabRcvMdl, 'enable', 'on');     
    %     set(handles.tabRfAffect, 'enable', 'on');        
    %     K = AMP .* exp(1i * PH);
    % else
    %     % ���� ������ ������ ��������� ��������, ����� ��������� ������� ���������
    %     set(handles.tabAcf, 'enable', 'off');
    %     set(handles.tabRcvParam, 'enable', 'off');
    %     set(handles.tabRcvMdl, 'enable', 'off'); 
    %     set(handles.tabRfAffect, 'enable', 'off');        
    %     K = [];
    % end
    
    % if (size(f, 2)>1)
    %     if (f0<min(f)) | (f0>max(f))  % ��� �� ������� ��������� ��������
    %         if get(handles.buttL1, 'value') == 1  % ������� �������� - L1
    %             buttL2_Callback(handles.buttL1, [], handles); % ������������� � L2
    %             %                msgbox(sprintf('��������������: ��������� ��������� �������������� �� ��������� � ��������� L1 (f0=%7.2f ���)', f0/1e6));
                
    %         else                                  % ������� �������� - L2
    %             buttL1_Callback(handles.buttL1, [], handles); % ������������� � L1
    %             %                msgbox(sprintf('��������������: ��������� ��������� �������������� �� ��������� � ��������� L2 (f0=%7.2f ���)', f0/1e6));
    %         end
    %     end
    % end



function edAcfPoints_Callback(hObject, eventdata, handles)
% hObject    handle to edAcfPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edAcfPoints as text
%        str2double(get(hObject,'String')) returns contents of edAcfPoints as a double

    s = get(hObject, 'string');
    n = str2num(s);
    if isempty(n)
        msgbox('� ������ ���� ���������� ������ ����� ������������� �����');
        set(hObject, 'string', 100);
    else
        if (n<1) | (n>1000000)
            msgbox('���������� ����� ������� ������ ������ � ��������� 1-1000000');
            set(hObject, 'string', 100);
        else
            set(hObject, 'string', round(n));
        end
    end

    

% --- Executes during object creation, after setting all properties.
function edAcfPoints_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edAcfPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edAcfTauMax_Callback(hObject, eventdata, handles)
% hObject    handle to edAcfTauMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edAcfTauMax as text
%        str2double(get(hObject,'String')) returns contents of edAcfTauMax as a double

    s = get(hObject, 'string');
    n = str2num(s);
    if isempty(n)
        msgbox('� ������ ���� ���������� ������ ����� � ������������');
        set(hObject, 'string', 100);
    else
        if (n<0.01) | (n>10000)
            msgbox('������� ������� ������ ������ � ��������� 0.01-10000 ��');
            set(hObject, 'string', 3000);
        else
            set(hObject, 'string', real(n));
        end
    end
    
    s = get(hObject, 'string');
    n = str2num(s);
    
    if get(handles.buttPT, 'value') == 1
        if (n < 2000) | (n>4000)
            msgbox('���� �������� ������� �� ����� ���������� � ������� ��������, ������ �������� ��������� ��������\n(�������� 2000-4000 ��)');
        end
    else
        if (n < 200) | (n>400)
            msgbox('���� �������� ������� �� ����� ���������� � ������� ��������, ������ �������� ��������� ��������\n (������� 200-400 ��)');
        end
    end
    


% --- Executes during object creation, after setting all properties.
function edAcfTauMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edAcfTauMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in buttAcfCalc.
function buttAcfCalc_Callback(hObject, eventdata, handles)
% hObject    handle to buttAcfCalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    global f K f0lit
    
    % if get(handles.buttL1, 'value') == 1
    f0 = 1602e6;
        if get(handles.buttPT, 'value') == 1
            sgn = 0;
        else
            sgn = 1;
        end
    % else
    %     %        f0 = 1247e6; 
    %     if get(handles.buttPT, 'value') == 1
    %         sgn = 2;
    %     else
    %         sgn = 3;
    %     end
    % end
    
    N =  str2num(get(handles.edAcfPoints, 'string'));          % ���������� �����
    TauMax = str2num(get(handles.edAcfTauMax, 'string'))*1e-9; % �������� �������, �����������
    
    % ����������, ������ ���
    %    H = waitbar(0, 'Please wait...');
    
    global hdl dt R R1
    hdl = handles;
    
    set(handles.panAcfProgr, 'visible', 'on');
    set(handles.lblAcfProgrDscr, 'visible', 'on');
    
    if get(handles.chkIFAoff, 'value') == 1
        IFAband = str2num(get(handles.edRcvIntBand, 'string'))*1e6;
    else
        IFAband = 0;
    end

    
    set(handles.buttSaveAcf, 'enable', 'off');

    if mod(N, 2) == 0
        N = N + 1;  % ����� ����� �������� �����, ����� dt=0 ���� ����������
    end
    
    [dt R R1] = AcfAnalyse(f, K, sgn, N, TauMax, f0lit, @AcfProgr, IFAband);
    set(handles.panAcfProgr, 'visible', 'off');
    set(handles.buttSaveAcf, 'enable', 'on');
    
       
    % ����� ���������� �� ������
    plot(handles.pltAcf, dt/1e-9, R/max(R), 'b', dt/1e-9, R1/max(R), 'r');
    grid(handles.pltAcf, 'on');
    
    % ������ ��������
    TauEst = sum(R1.*dt) / sum(R1);
    TauEst*1e9
    
    
    set(handles.lblTauEst, 'String', sprintf('%10.3f', abs(TauEst*1e9)));
    
    set(handles.pltAcf, 'ButtonDownFcn',  @pltAcf_ButtonDownFcn);
    
function AcfProgr(x, dscr)
    global hdl
    
    set(hdl.lblAcfProgr, 'Position', [0 0.5 min(0.01+x, 1) 0.45]);
    set(hdl.lblAcfProgrDscr, 'string', dscr);
    drawnow


% --- Executes on selection change in popRcvType.
function popRcvType_Callback(hObject, eventdata, handles)
% hObject    handle to popRcvType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popRcvType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popRcvType

    global UsrRcvParams PredRcv

    vr = get(hObject, 'value');
    
    set(handles.edRcvIntBand, 'string', UsrRcvParams(vr, 1));
    set(handles.edRcvPllBand, 'string', UsrRcvParams(vr, 2));
    set(handles.edRcvDllBand, 'string', UsrRcvParams(vr, 3));

    if vr < PredRcv+1 % ��������� ������ PredRcv ��������� ������ � ���������, ������ �� ������
        set(handles.edRcvIntBand, 'enable', 'off');
        set(handles.edRcvPllBand, 'enable', 'off');        
        set(handles.edRcvDllBand, 'enable', 'off');                
    else    
        if get(handles.chkIFAoff, 'value') == 1
            set(handles.edRcvIntBand, 'enable', 'on');
        else
            set(handles.edRcvIntBand, 'enable', 'off');
        end
        set(handles.edRcvPllBand, 'enable', 'on');        
        set(handles.edRcvDllBand, 'enable', 'on');                
    end


% --- Executes during object creation, after setting all properties.
function popRcvType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popRcvType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edRcvIntBand_Callback(hObject, eventdata, handles)
% hObject    handle to edRcvIntBand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edRcvIntBand as text
%        str2double(get(hObject,'String')) returns contents of edRcvIntBand as a double

global UsrRcvParams

vr = get(handles.popRcvType, 'value');
if vr > 2
    s = get(hObject, 'string');
    n = str2num(s);
    if isempty(n)
        msgbox('� ������ ���� ���������� ������ ������������� �����');
        set(hObject, 'string', 10);
    else
        if (n<0.5) | (n>50)
            msgbox('������ ��� ������ ������ � ��������� 0.5-50 ���');
            set(hObject, 'string', 10);
        else
            set(hObject, 'string', real(n));
        end
        UsrRcvParams(vr, 1) = real(n);
    end
end

% --- Executes during object creation, after setting all properties.
function edRcvIntBand_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edRcvIntBand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edRcvPllBand_Callback(hObject, eventdata, handles)
% hObject    handle to edRcvPllBand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edRcvPllBand as text
%        str2double(get(hObject,'String')) returns contents of edRcvPllBand as a double

global UsrRcvParams

vr = get(handles.popRcvType, 'value');
if vr > 2
    s = get(hObject, 'string');
    n = str2num(s);
    if isempty(n)
        msgbox('� ������ ���� ���������� ������ ������������� �����');
        set(hObject, 'string', 30);
    else
        if (n<1) | (n>50)
            msgbox('������ ��� ������ ������ � ��������� 1-50 ��');
            set(hObject, 'string', 10);
        else
            set(hObject, 'string', real(n));
        end
        UsrRcvParams(vr, 2) = real(n);
    end
end


% --- Executes during object creation, after setting all properties.
function edRcvPllBand_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edRcvPllBand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edRcvDllBand_Callback(hObject, eventdata, handles)
% hObject    handle to edRcvDllBand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edRcvDllBand as text
%        str2double(get(hObject,'String')) returns contents of edRcvDllBand as a double
global UsrRcvParams

vr = get(handles.popRcvType, 'value');
if vr > 2
    s = get(hObject, 'string');
    n = str2num(s);
    if isempty(n)
        msgbox('� ������ ���� ���������� ������ ������������� �����');
        set(hObject, 'string', 3);
    else
        if (n<0.1) | (n>20)
            msgbox('������ ��� ������ ������ � ��������� 0.1-20 ��');
            set(hObject, 'string', 10);
        else
            set(hObject, 'string', real(n));
        end
        UsrRcvParams(vr, 3) = real(n);
    end
end


% --- Executes during object creation, after setting all properties.
function edRcvDllBand_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edRcvDllBand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tabRcvMdl.
function tabRcvMdl_Callback(hObject, eventdata, handles)
% hObject    handle to tabRcvMdl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    set(handles.tabSrcData, 'value', 0);
    set(handles.tabAcf, 'value', 0);
    set(handles.tabRcvParam, 'value', 0);
    set(handles.tabRcvMdl, 'value', 1);
    set(handles.tabRfAffect, 'value', 0);

    set(handles.panLoadSrcData, 'visible', 'off');
    set(handles.panAcf, 'visible', 'off');
    set(handles.panRcvParam, 'visible', 'off');
    set(handles.panRcvMdl, 'visible', 'on');
    set(handles.panRfAffect, 'visible', 'off');        


% --- Executes on button press in buttRcvMdlStart.
function buttRcvMdlStart_Callback(hObject, eventdata, handles)
% hObject    handle to buttRcvMdlStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global f K sgn f0lit


% if get(handles.buttL1, 'value') == 1
f0 = 1602e6;
    if get(handles.buttPT, 'value') == 1
        sgn = 0;
    else
        sgn = 1;
    end
% else
%     %    f0 = 1247e6; 
%     if get(handles.buttPT, 'value') == 1
%         sgn = 2;
%     else
%         sgn = 3;
%     end
% end

s = get(handles.edRcvMdlTime, 'string');
TMax = str2num(s);

s = get(handles.edRcvPllBand, 'string');
dff = str2num(s);

s = get(handles.edRcvDllBand, 'string');
dft = str2num(s);

global hdl RcvMdlStopFlag c
hdl = handles;
        
RcvMdlDisable();

set(handles.edRcvMdlMFreq, 'string', 0);
set(handles.edRcvMdlDFreq, 'string', 0);
set(handles.edRcvMdlMTau, 'string', 0);
set(handles.edRcvMdlDTau, 'string', 0);

set(handles.lblErrMsg, 'String', '');

RcvMdlStopFlag = 0;
set(handles.panRcvMdlProgr, 'visible', 'on');
set(handles.lblRcvMdlProgrDscr, 'visible', 'on');

if get(handles.chkIFAoff, 'value') == 1
    IFAband = str2num(get(handles.edRcvIntBand, 'string'))*1e6;
else
    IFAband = 0;
end

q = 10^(str2num(get(handles.edRcvQ, 'string'))/10);
V0 = str2num(get(handles.edRcvVel, 'string'));
A0 = str2num(get(handles.edRcvAcc, 'string'));

if get(handles.buttTauMode, 'value')==1
    TrkMode = 0;
else
    TrkMode = 1;
end

global dt1 ErrF ErrT

set(handles.buttSaveErrFreq, 'enable', 'off');
set(handles.buttSaveErrTau, 'enable', 'off');


[dt1 ErrP ErrF ErrT] = TrackModel(f, K, q, sgn, TMax, f0lit, dff, dft, @RcvMdlProgr, handles.pltErrFreq, handles.pltErrTau, IFAband, V0, A0, handles.lblErrMsg, TrkMode);
set(handles.panRcvMdlProgr, 'visible', 'off');

set(handles.buttSaveErrFreq, 'enable', 'on');
set(handles.buttSaveErrTau, 'enable', 'on');


set(handles.pltErrFreq, 'ButtonDownFcn', @pltErrFreq_ButtonDownFcn);
set(handles.pltErrTau, 'ButtonDownFcn',  @pltErrTau_ButtonDownFcn);

RcvMdlEnable();

if size(ErrF, 2) < 1000
    set(handles.lblErrMsg, 'string', '��� ����� ���������� - �� ������� ���������� �������');
    set(handles.edRcvMdlMFreq, 'string', mean(ErrF(:))/2/pi/f0lit*c);
    set(handles.edRcvMdlDFreq, 'string', std(ErrF(:))/2/pi/f0lit*c);
    set(handles.edRcvMdlMTau, 'string', mean(ErrT(:)));
    set(handles.edRcvMdlDTau, 'string', std(ErrT(:)));
    
else
    i = 100:size(ErrF, 2);
    set(handles.edRcvMdlMFreq, 'string', mean(ErrF(i))/2/pi/f0lit*c);
    set(handles.edRcvMdlDFreq, 'string', std(ErrF(i))/2/pi/f0lit*c);
    i = 1000:size(ErrF, 2);
    set(handles.edRcvMdlMTau, 'string', mean(ErrT(i)));
    set(handles.edRcvMdlDTau, 'string', std(ErrT(i)));
end
    
    
function RcvMdlProgr(x, dscr)
    global hdl
    
    set(hdl.lblRcvMdlProgr, 'Position', [0 0.5 min(0.01+x, 1) 0.45]);
    set(hdl.lblRcvMdlProgrDscr, 'string', dscr);
    drawnow



function edRcvMdlTime_Callback(hObject, eventdata, handles)
% hObject    handle to edRcvMdlTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edRcvMdlTime as text
%        str2double(get(hObject,'String')) returns contents of edRcvMdlTime as a double

    s = get(hObject, 'string');
    n = str2num(s);
    if isempty(n)
        msgbox('� ������ ���� ���������� ������ ����� ������������� �����');
        set(hObject, 'string', 1);
    else
        if (n<0.01) | (n>10000)
            msgbox('����� ������������� ������ ���� � ��������� 0.01-10000 �');
            set(hObject, 'string', 1);
        else
            set(hObject, 'string', real(n));
        end
    end


% --- Executes during object creation, after setting all properties.
function edRcvMdlTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edRcvMdlTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edRcvMdlMFreq_Callback(hObject, eventdata, handles)
% hObject    handle to edRcvMdlMFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edRcvMdlMFreq as text
%        str2double(get(hObject,'String')) returns contents of edRcvMdlMFreq as a double


% --- Executes during object creation, after setting all properties.
function edRcvMdlMFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edRcvMdlMFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edRcvMdlDFreq_Callback(hObject, eventdata, handles)
% hObject    handle to edRcvMdlDFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edRcvMdlDFreq as text
%        str2double(get(hObject,'String')) returns contents of edRcvMdlDFreq as a double


% --- Executes during object creation, after setting all properties.
function edRcvMdlDFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edRcvMdlDFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edRcvMdlMTau_Callback(hObject, eventdata, handles)
% hObject    handle to edRcvMdlMTau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edRcvMdlMTau as text
%        str2double(get(hObject,'String')) returns contents of edRcvMdlMTau as a double


% --- Executes during object creation, after setting all properties.
function edRcvMdlMTau_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edRcvMdlMTau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edRcvMdlDTau_Callback(hObject, eventdata, handles)
% hObject    handle to edRcvMdlDTau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edRcvMdlDTau as text
%        str2double(get(hObject,'String')) returns contents of edRcvMdlDTau as a double


% --- Executes during object creation, after setting all properties.
function edRcvMdlDTau_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edRcvMdlDTau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edRcvMdlTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edRcvMdlTime as text
%        str2double(get(hObject,'String')) returns contents of edRcvMdlTime as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edRcvMdlTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in buttRcvMdlStop.
function buttRcvMdlStop_Callback(hObject, eventdata, handles)
% hObject    handle to buttRcvMdlStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RcvMdlStopFlag
RcvMdlStopFlag = 1;

set(handles.buttRcvMdlStart, 'Enable', 'on');
set(handles.buttRcvMdlStop, 'Enable', 'off');

RcvMdlEnable();


function RcvMdlDisable()
    global RcvMdlDisables hdl

    set(hdl.buttRcvMdlStart, 'Enable', 'off');
    set(hdl.buttRcvMdlStop, 'Enable', 'on');
    set(hdl.chkIFAoff, 'Enable', 'off');
    
    for i=1:size(RcvMdlDisables, 1)
        set(RcvMdlDisables(i),     'Enable', 'off');
    end
    
    
function RcvMdlEnable()
    global RcvMdlDisables hdl PredRcv

    set(hdl.buttRcvMdlStart, 'Enable', 'on');
    set(hdl.buttRcvMdlStop, 'Enable', 'off');
    set(hdl.chkIFAoff, 'Enable', 'off');
    
    for i=1:size(RcvMdlDisables, 1)
        set(RcvMdlDisables(i),     'Enable', 'on');
    end
    
    popRcvType_Callback(hdl.popRcvType, [], hdl);
    


% --- Executes on button press in chkIFAoff.
function chkIFAoff_Callback(hObject, eventdata, handles)
% hObject    handle to chkIFAoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkIFAoff

    global PredRcv
    
    if get(hObject, 'Value') == 1
        set(hObject, 'String', '���');
        if get(handles.popRcvType, 'value') < PredRcv
            set(handles.edRcvIntBand, 'Enable', 'off');
        else
            set(handles.edRcvIntBand, 'Enable', 'on');
        end
    else
        set(hObject, 'String', '����');
        set(handles.edRcvIntBand, 'Enable', 'off');
    end
    



function edRcvQ_Callback(hObject, eventdata, handles)
% hObject    handle to edRcvQ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edRcvQ as text
%        str2double(get(hObject,'String')) returns contents of edRcvQ as a double

    s = get(hObject, 'string');
    n = str2num(s);
    if isempty(n)
        msgbox('� ������ ���� ���������� ������ ������������� �����');
        set(hObject, 'string', 45);
    else
        if (n<25) | (n>60)
            msgbox('��������� �������� ������� � ������������ ��������� ���� ������ ������ � �������� 25-60 ����');
            set(hObject, 'string', 45);
        else
            set(hObject, 'string', real(n));
        end
    end


% --- Executes during object creation, after setting all properties.
function edRcvQ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edRcvQ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edRcvVel_Callback(hObject, eventdata, handles)
% hObject    handle to edRcvVel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edRcvVel as text
%        str2double(get(hObject,'String')) returns contents of edRcvVel as a double


    s = get(hObject, 'string');
    n = str2num(s);
    if isempty(n)
        msgbox('� ������ ���� ���������� ������ �����');
        set(hObject, 'string', 0);
    else
        if (n<-20000) | (n>20000)
            msgbox('�������� ������ ������ � �������� +- 20000 �/� ');
            set(hObject, 'string', 0);
        else
            set(hObject, 'string', real(n));
        end
    end


% --- Executes during object creation, after setting all properties.
function edRcvVel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edRcvVel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edRcvAcc_Callback(hObject, eventdata, handles)
% hObject    handle to edRcvAcc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edRcvAcc as text
%        str2double(get(hObject,'String')) returns contents of edRcvAcc as a double
    s = get(hObject, 'string');
    n = str2num(s);
    if isempty(n)
        msgbox('� ������ ���� ���������� ������ �����');
        set(hObject, 'string', 0);
    else
        if (n<-200) | (n>200)
            msgbox('��������� ������ ������ � �������� +- 200 �/�^2 ');
            set(hObject, 'string', 0);
        else
            set(hObject, 'string', real(n));
        end
    end


% --- Executes during object creation, after setting all properties.
function edRcvAcc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edRcvAcc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popLit.
function popLit_Callback(hObject, eventdata, handles)
% hObject    handle to popLit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popLit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popLit

global f0lit

ind = get(handles.popLit, 'Value');
Lits = get(handles.popLit, 'String');
%if get(handles.buttL1, 'Value') == 1 % L1
%    f0lit = 1602e6 + 562.5e3 * str2num(Lits{ind});
   f0lit = 1e6*str2num(Lits{ind});
% else  % L2
%     f0lit = 1246e6 + 437.5e3 * str2num(Lits{ind});
% end
set(handles.lblLit, 'String', sprintf('f0=%7.2f ��', f0lit/1e6));


% --- Executes during object creation, after setting all properties.
function popLit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popLit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on mouse press over axes background.
function pltErrFreq_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pltErrFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
axesfig = figure;
% Copy the axes and size it to the figure
axescopy = copyobj(hObject,axesfig);
set(axescopy,'Units','Normalized',...
              'Position',[.05,.10,.90,.8])

% --- Executes on mouse press over axes background.
function pltErrTau_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pltErrFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
axesfig = figure;
% Copy the axes and size it to the figure
axescopy = copyobj(hObject,axesfig);
set(axescopy,'Units','Normalized',...
              'Position',[.05,.10,.90,.8])


% --- Executes on mouse press over axes background.
function pltAcf_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pltErrFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
axesfig = figure;
% Copy the axes and size it to the figure
axescopy = copyobj(hObject,axesfig);
set(axescopy,'Units','Normalized',...
              'Position',[.05,.10,.90,.8])


% --- Executes on mouse press over axes background.
function pltAFR_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pltErrFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
axesfig = figure;
% Copy the axes and size it to the figure
axescopy = copyobj(hObject,axesfig);
set(axescopy,'Units','Normalized',...
              'Position',[.05,.10,.90,.8])


% --- Executes on mouse press over axes background.
function pltPFR_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pltErrFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
axesfig = figure;
% Copy the axes and size it to the figure
axescopy = copyobj(hObject,axesfig);
set(axescopy,'Units','Normalized',...
              'Position',[.05,.10,.90,.8])



% --- Executes on button press in buttTauMode.
function buttTauMode_Callback(hObject, eventdata, handles)
% hObject    handle to buttTauMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of buttTauMode


    set(handles.buttPhaseMode, 'value', 0);
    
% --- Executes on button press in buttPhaseMode.
function buttPhaseMode_Callback(hObject, eventdata, handles)
% hObject    handle to buttPhaseMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of buttPhaseMode

    set(handles.buttTauMode, 'value', 0);


% --- Executes on button press in buttSaveAcf.
function buttSaveAcf_Callback(hObject, eventdata, handles)
% hObject    handle to buttSaveAcf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    global dt R R1
    [sFileName sPathName] = uiputfile('*.csv');    

    fid = fopen([sPathName sFileName], 'w');    
    fprintf(fid, '%g,%g\n', [dt; R1]);
    fclose(fid);
    
    sFileName2 = [sFileName(1:size(sFileName, 2)-4) '_ideal' sFileName(size(sFileName, 2)-3:size(sFileName, 2))];
    fid = fopen([sPathName sFileName2], 'w');    
    fprintf(fid, '%g,%g\n', [dt; R]);
    fclose(fid);
    
    
    

% --- Executes on button press in buttSaveErrFreq.
function buttSaveErrFreq_Callback(hObject, eventdata, handles)
% hObject    handle to buttSaveErrFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    global dt1 ErrF
    [sFileName sPathName] = uiputfile('*.csv');    

    fid = fopen([sPathName sFileName], 'w');    
    fprintf(fid, '%g,%g\n', [dt1; ErrF]);
    fclose(fid);
    

    

% --- Executes on button press in buttSaveErrTau.
function buttSaveErrTau_Callback(hObject, eventdata, handles)
% hObject    handle to buttSaveErrTau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    global dt1 ErrT
    [sFileName sPathName] = uiputfile('*.csv');    

    fid = fopen([sPathName sFileName], 'w');    
    fprintf(fid, '%g,%g\n', [dt1; ErrT]);
    fclose(fid);



function edRfAffectTime_Callback(hObject, eventdata, handles)
% hObject    handle to edRfAffectTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edRfAffectTime as text
%        str2double(get(hObject,'String')) returns contents of edRfAffectTime as a double
    s = get(hObject, 'string');
    n = str2num(s);
    if isempty(n)
        msgbox('� ������ ���� ���������� ������ ����� ������������� �����');
        set(hObject, 'string', 1000);
    else
        if (n<1) | (n>1000)
            msgbox('���������� ����� ������� ������ ������ � ��������� 1-1000');
            set(hObject, 'string', 1000);
        else
            set(hObject, 'string', round(n));
        end
    end


% --- Executes during object creation, after setting all properties.
function edRfAffectTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edRfAffectTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in buttRfAffectStart.
function buttRfAffectStart_Callback(hObject, eventdata, handles)
% hObject    handle to buttRfAffectStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global f K f0lit
    
    % if get(handles.buttL1, 'value') == 1
    f0 = 1602e6;
        if get(handles.buttPT, 'value') == 1
            sgn = 0;
        else
            sgn = 1;
        end
    % else
    %     %        f0 = 1247e6; 
    %     if get(handles.buttPT, 'value') == 1
    %         sgn = 2;
    %     else
    %         sgn = 3;
    %     end
    % end
    
    TauMax =  1e-6*str2num(get(handles.edRfAffectTime, 'string'));          % ���������� �����
                                                                   %    TauMax = str2num(get(handles.edAcfTauMax, 'string'))*1e-9; % �������� �������, �����������
    N = 1000;
    
    % ����������, ������ ���
    %    H = waitbar(0, 'Please wait...');
    
    global hdl dt R R1
    hdl = handles;
    
    set(handles.panAcfProgr, 'visible', 'on');
    set(handles.lblAcfProgrDscr, 'visible', 'on');
    
    if get(handles.chkIFAoff, 'value') == 1
        IFAband = str2num(get(handles.edRcvIntBand, 'string'))*1e6;
    else
        IFAband = 0;
    end

    
    %    set(handles.buttSaveAcf, 'enable', 'off');
    
    [t s y PSP] = RfAffectAnalyse(f, K, sgn, N, TauMax, f0lit, @AcfProgr, IFAband);
    set(handles.panAcfProgr, 'visible', 'off');
    set(handles.buttSaveAcf, 'enable', 'on');
    
       
    % ����� ���������� �� ������
    plot(handles.pltRfAffect, t/1e-9, real(s), 'b', t/1e-9, real(y)/std(real(y))/sqrt(2)+2.1, 'k', t/1e-9, PSP, 'r');
    grid(handles.pltRfAffect, 'on');
    
    % % ������ ��������
    % TauEst = sum(R1.*dt) / sum(R1.*R1);
    
    % set(handles.lblTauEst, 'String', sprintf('%7.2f', abs(TauEst*1e9)));
    
    set(handles.pltRfAffect, 'ButtonDownFcn',  @pltRfAffect_ButtonDownFcn);    

    

% --- Executes on button press in buttPTRf.
function buttPTRf_Callback(hObject, eventdata, handles)
% hObject    handle to buttPTRf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of buttPTRf
    set(handles.buttPT, 'Value', 1);
    set(handles.buttVT, 'Value', 0);
    set(handles.buttPTRf, 'Value', 1);
    set(handles.buttVTRf, 'Value', 0);


% --- Executes on button press in buttVTRf.
function buttVTRf_Callback(hObject, eventdata, handles)
% hObject    handle to buttVTRf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of buttVTRf
    set(handles.buttPT, 'Value', 0);
    set(handles.buttVT, 'Value', 1);
    set(handles.buttPTRf, 'Value', 0);
    set(handles.buttVTRf, 'Value', 1);


% --- Executes on selection change in popLitRf.
function popLitRf_Callback(hObject, eventdata, handles)
% hObject    handle to popLitRf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popLitRf contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popLitRf
% hObject    handle to popLit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popLit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popLit

global f0lit

ind = get(handles.popLitRf, 'Value');
Lits = get(handles.popLitRf, 'String');
%if get(handles.buttL1, 'Value') == 1 % L1
%    f0lit = 1602e6 + 562.5e3 * str2num(Lits{ind});
    f0lit = 1e6 * str2num(Lits{ind});
% else  % L2
%     f0lit = 1246e6 + 437.5e3 * str2num(Lits{ind});
% end
set(handles.lblLit, 'String', sprintf('f0=%7.2f ��', f0lit/1e6));  % !!! lblLit !!!



% --- Executes during object creation, after setting all properties.
function popLitRf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popLitRf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pltRfAffect_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pltErrFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    
    456
    
axesfig = figure;
% Copy the axes and size it to the figure
axescopy = copyobj(hObject,axesfig);
set(axescopy,'Units','Normalized',...
              'Position',[.05,.10,.90,.8])


% --- Executes on button press in tabRfAffect.
function tabRfAffect_Callback(hObject, eventdata, handles)
% hObject    handle to tabRfAffect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tabRfAffect
    set(handles.tabSrcData, 'value', 0);
    set(handles.tabAcf, 'value', 0);
    set(handles.tabRcvParam, 'value', 0);
    set(handles.tabRcvMdl, 'value', 0);
    set(handles.tabRfAffect, 'value', 1);

    set(handles.panLoadSrcData, 'visible', 'off');
    set(handles.panAcf, 'visible', 'off');
    set(handles.panRcvParam, 'visible', 'off');
    set(handles.panRcvMdl, 'visible', 'off');
    set(handles.panRfAffect, 'visible', 'on');    
