function varargout = playing_with_OFDM(varargin)
% PLAYING_WITH_OFDM MATLAB code for playing_with_OFDM.fig
%      PLAYING_WITH_OFDM, by itself, creates a new PLAYING_WITH_OFDM or raises the existing
%      singleton*.
%
%      H = PLAYING_WITH_OFDM returns the handle to a new PLAYING_WITH_OFDM or the handle to
%      the existing singleton*.
%
%      PLAYING_WITH_OFDM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLAYING_WITH_OFDM.M with the given input arguments.
%
%      PLAYING_WITH_OFDM('Property','Value',...) creates a new PLAYING_WITH_OFDM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before playing_with_OFDM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to playing_with_OFDM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help playing_with_OFDM

% Last Modified by GUIDE v2.5 12-Oct-2016 17:01:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @playing_with_OFDM_OpeningFcn, ...
    'gui_OutputFcn',  @playing_with_OFDM_OutputFcn, ...
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


% --- Executes just before playing_with_OFDM is made visible.
function playing_with_OFDM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to playing_with_OFDM (see VARARGIN)

% Choose default command line output for playing_with_OFDM
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes playing_with_OFDM wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = playing_with_OFDM_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double

% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
global snrValNew;
try
    snrValNew = eval(get(hObject,'String'));
catch
    set(handles.text18, 'String','enter a valid SNR');
end

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
global pathGainsNew;
try
    pathGainsNew = eval(get(hObject,'String'));
catch
    set(handles.text18, 'String','enter valid path gains');
end

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
global fdopNew;
try
    fdopNew = eval(get(hObject,'String'));
catch
    set(handles.text18, 'String','enter a valid doppler frequency');
end

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
global fsNew;
try
    fsNew = eval(get(hObject,'String')); % sampling frequency
catch
    set(handles.text18, 'String','enter a valid sampling frequency');
end

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
global bwSigNew;
try
    bwSigNew = eval(get(hObject,'String')); % signal bandwidth
catch
    set(handles.text18, 'String','enter a valid bandwidth');
end

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
global fsubNew;
try
    fsubNew = eval(get(hObject,'String')); % sub-carrier spacing
catch
    set(handles.text18, 'String','enter a valid sub-carrier spacing');
end

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
global symRateNew; % symbol rate
try
    symRateNew = eval(get(hObject,'String')); %  symbol rate
catch
    set(handles.text18, 'String','enter a valid OFDM symbol rate');
end

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
global NcpNew; % number of CP samples
global fsubNew;
global bwSigNew;
Nsc = round(bwSigNew/fsubNew);
try
    NcpTmp = eval(get(hObject,'String')); %  number of CP samples at sampling rate Fs
catch
    set(handles.text18, 'String','enter a valid CP length');
end
if NcpTmp >= Nsc
    set(handles.text18, 'String','choose a smaller CP length');
else
    NcpNew = NcpTmp;
    set(handles.text18, 'String','');
end

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
global modTypeStrNew;
contents = cellstr(get(hObject,'String'));
modTypeStrNew = contents{get(hObject,'Value')};

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double
global dpSpeed;
dpSpeed = eval(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double
global numOfdmSymsOnGridNew;
try
    numOfdmSymsOnGridNew = eval(get(hObject,'String'));
catch
    set(handles.text18, 'String','enter a valid number of OFDM symbols to display');
end

% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dbstop if error;

global setDefaultFlag;

if isempty(setDefaultFlag)
    set_default_vals();
end

global bwSig;
global fs;
global fsChan;
global fsub;
global symRate;
global fdop;
global modTypeStr;
global pathGains;
global pathDelays;
global Ncp;
global numOfdmSymsOnGrid;
global viewTimeDomainWaveform;
global viewSubCarriers;
global pauseSim;
global exitButton;
global viewChEst3D;
global viewChEstMagPhase;
global useIdealChEstForEq;
global snrVal;
global cfoVal;
global viewEVM;


%%
chStruct.doppFilCoeff = get_dopp_filter(fdop, fsChan);
chanState = [];
sigState = [];
awgnDoppLastSamp = [];
fdopPrev = fdop;
viewSubCarrierOsr = 10;
numViewSubCarriers = 10;
fsubSymRatioPrev = fsub/symRate;
firstTimePlotFlg = 0;
fsPrev = fs;
NcpPrev = Ncp;
bwSigPrev = bwSig;
fsubPrev = fsub;

%doppFilOpArr = [];
pilotSpacing = 3;
viewTimeDomainWaveformPrev = 0;
viewSubCarriersPrev = 0;
viewChEst3DPrev = 0;
viewChEstMagPhasePrev = 0;

% generate FFT input based on modulation type
for loopCnt = 1:10000
    try
        if loopCnt == 1 || ~isequal(fsPrev, fs) || ~isequal(NcpPrev, Ncp) || ~isequal(fsubPrev, fsub) || ~isequal(bwSigPrev, bwSig)
            Nsc = round(bwSig/fsub); % number of sub-carriers = bandwidth/sub-carrier spacing
            Nsc = double(mod(Nsc, 2) == 0) * (Nsc-1) + double(mod(Nsc, 2) == 1) * Nsc;
            fftOpArr = zeros(Nsc, numOfdmSymsOnGrid); % raw received OFDM symbols (Resource Elements)
            fftOpEqArr = zeros(Nsc, numOfdmSymsOnGrid); % equalized OFDM symbols (Zero-forcing)
            fftOpEqVec = zeros(1, Nsc*numOfdmSymsOnGrid);
            chEstInterpArr = zeros(Nsc, numOfdmSymsOnGrid);
            txSigArr = zeros(1, (Nsc+Ncp) * numOfdmSymsOnGrid);
            idealChanFreqDomainArr = zeros(Nsc, numOfdmSymsOnGrid); % ideal channel for "numOfdmSymsOnGrid" OFDM symbols
            fsPrev = fs;
            NcpPrev = Ncp;
            fsubPrev = fsub;
            bwSigPrev = bwSig;
            cfoValPrev = 0;
            fsChan = ceil(fs/(Nsc+Ncp));
        end
        if pauseSim == 1
            set(handles.pushbutton5,'string','RESUME','enable','on');
            if exitButton == 1;
                close(handles.figure1);
                return;
            end
            pause(1);continue;
        else
            set(handles.pushbutton5,'string','PAUSE','enable','on');
        end
        
        switch modTypeStr
            case 'BPSK'
                fftIp = (2*round(rand(1, Nsc))-1);
                evmConstRef = unique((2*round(rand(1, 1e2))-1));
                txPwr = 1;
            case 'QPSK'
                fftIp = (2*round(rand(1, Nsc))-1) + 1j*(2*round(rand(1, Nsc)) -1);
                evmConstRef = unique((2*round(rand(1, 1e2))-1) + 1j*(2*round(rand(1, 1e2))-1));
                txPwr = 2;
            case '16QAM'
                fftIp = (2*round(3*rand(1, Nsc))-3) + 1j*(2*round(3*rand(1, Nsc)) -3);
                evmConstRef = unique((2*round(3*rand(1, 1e2))-3) + 1j*(2*round(3*rand(1, 1e2))-3));
                txPwr = 10;
            case '64QAM'
                fftIp = (2*round(7*rand(1, Nsc))-7) + 1j*(2*round(7*rand(1, Nsc)) -7);
                evmConstRef = unique((2*round(7*rand(1, 1e3))-7) + 1j*(2*round(7*rand(1, 1e3))-7));
                txPwr = 42;
            case '256QAM'
                fftIp = (2*round(15*rand(1, Nsc))-15) + 1j*(2*round(15*rand(1, Nsc)) -15);
                evmConstRef = unique((2*round(15*rand(1, 1e3))-15) + 1j*(2*round(15*rand(1, 1e3))-15));
                txPwr = 170;
            otherwise % DEFAULT is QPSK
                fftIp = (2*round(rand(1, Nsc))-1) + 1j*(2*round(rand(1, Nsc)) -1);
                evmConstRef = unique((2*round(rand(1, 1e2))-1) + 1j*(2*round(rand(1, 1e2))-1));
                txPwr = 2;
        end
        
        % insert pilot symbols
        %chEstIdeal = (1+1j)/sqrt(2)*sqrt(txPwr);
        chEstIdeal = max(real(evmConstRef)) + 1j*max(imag(evmConstRef));
        fftIp([1:pilotSpacing:end]) = chEstIdeal;
        
        fftIpIntf = fftIp*sinc(1/symRate*([0:fsub:fsub*Nsc-fsub]'*ones(1,Nsc)-ones(Nsc, 1)*[0:1:Nsc-1]*fsub)).';
        
        if fdop ~= fdopPrev
            chStruct.doppFilCoeff = get_dopp_filter(fdop, fsChan);
        end
        fdopPrev = fdop;
        
        % insert CP and take IFFT
        txSig = ifft([fftIpIntf])*sqrt(Nsc);%/Nsc;
        txSig = [txSig(end-Ncp+1:end) txSig];
        
        txSigArr(1:(Ncp+Nsc)*(numOfdmSymsOnGrid-1)) = txSigArr(Ncp+Nsc+1:end);
        txSigArr((Ncp+Nsc)*(numOfdmSymsOnGrid-1)+1:end) = txSig;
        
        % Apply multipath channel
        [multipathChanOp, sigState, chanState, idealChan, doppFilResampled, awgnDoppLastSamp] = myrayleighchan(txSig, chStruct, pathGains, pathDelays, fs, fsChan, sigState, chanState, awgnDoppLastSamp);
        
        % Apply noise
        chPwr = sum(10.^(pathGains/10));
        noiseVec = sqrt(chPwr*txPwr)*(1/sqrt(2))*(randn(size(multipathChanOp)) + 1j * randn(size(multipathChanOp)))*10^(-snrVal/20);
        
        rxSigNoise = multipathChanOp + noiseVec;
        
        % Apply CFO and sample offset
        rxSigNoiseCfo = rxSigNoise.*exp(1j*2*pi*[0:1:Nsc+Ncp-1]*cfoVal/fs + cfoValPrev);
        cfoValPrev = 1j*2*pi*(Nsc+Ncp-1)*cfoVal/fs + cfoValPrev;
        
        % CP removal
        % FFT
        fftOp = fft(rxSigNoiseCfo(Ncp+1:end));
        chEstRaw = fftOp(1:pilotSpacing:end)/chEstIdeal;
        pilotIndices = [1:pilotSpacing:ceil(Nsc/pilotSpacing)*pilotSpacing];
        try
            chEstInterp = interp1(pilotIndices,chEstRaw, [1:1:max(pilotIndices)]);
        catch
            chEstInterp = chEstRaw;
        end
        chEstInterp(max(pilotIndices)+1:Nsc) = chEstInterp(max(pilotIndices));
        
        % ideal channel frequency domain
        idealChanFreqDomain = fft(idealChan(1:end-Ncp))*sqrt(Nsc);
        
        % equalization
        if useIdealChEstForEq == 1
            fftOpEq = fftOp./idealChanFreqDomain;
        else
            fftOpEq = fftOp./chEstInterp;
        end
        
        idealChanFreqDomainArr(:, 1:end-1) = idealChanFreqDomainArr(:, 2:end);
        idealChanFreqDomainArr(:, end) = idealChanFreqDomain.';
        
        fftOpArr(:, 1:end-1) = fftOpArr(:, 2:end);
        fftOpArr(:, end) = fftOp.';
        fftOpEqArr(:, 1:end-1) = fftOpEqArr(:, 2:end);
        fftOpEqArr(:, end) = fftOpEq.';
        fftOpEqVec(1:(numOfdmSymsOnGrid-1)*Nsc) = fftOpEqVec([1:(numOfdmSymsOnGrid-1)*Nsc] + Nsc);
        fftOpEqVec([1:Nsc]+(numOfdmSymsOnGrid-1)*Nsc) = fftOpEq;
        chEstInterpArr(:, 1:end-1) = chEstInterpArr(:, 2:end);
        chEstInterpArr(:, end) = chEstInterp.';
        
        
        surf(handles.axes1, abs(idealChanFreqDomainArr)); zLimit = max(max(abs(idealChanFreqDomainArr)));
        
        xlabel(handles.axes1, 'OFDM symbol');ylabel(handles.axes1, 'sub-carrier');zlabel(handles.axes1, 'Amplitude');
        set(handles.axes1, 'zlim', [0 zLimit]);
        %drawnow();
        plot(handles.axes2, fftOpEqVec, 'k*');
        xlabel(handles.axes2, 'Real');ylabel(handles.axes2, 'Imaginary');
        set(handles.axes2, 'xlim', [-16 16], 'ylim', [-16 16]);
        drawnow();
        
        if viewEVM == 1
            [evmVal] = evm_compute(fftOpEq, evmConstRef);
            set(handles.edit11,'string',evmVal);
        else
            set(handles.edit11,'string','NA');
        end
        
        % VIEW TIME DOMAIN WAVEFORM
        if viewTimeDomainWaveform == 1
            figure(1);plot([0:1/fs:(Nsc+Ncp)*numOfdmSymsOnGrid/fs-1/fs], abs(txSigArr), 'b');xlabel('time (seconds)');ylabel('Amplitude');title('OFDM time domain waveform (magnitude) at fs');
        else
            if viewTimeDomainWaveform ~= viewTimeDomainWaveformPrev
                close(figure(1));
            end
        end
        viewTimeDomainWaveformPrev = viewTimeDomainWaveform;
        
        % VIEW SUB-CARRIERS (WIDTH AND SPACING RELATION)
        if viewSubCarriers == 1
            fsubSymRatio = fsub/symRate;
            if fsubSymRatio ~= fsubSymRatioPrev
                firstTimePlotFlg = 0;
            end
            % the sub-carrier plot does not change for every iteration
            % it changes only when the symbol rate or the sub-carrier spacing is changed
            if firstTimePlotFlg == 0
                sincMat = sinc(1/symRate*([0:fsub/viewSubCarrierOsr:numViewSubCarriers*fsub-fsub/viewSubCarrierOsr]'*ones(1,numViewSubCarriers)-ones(numViewSubCarriers*viewSubCarrierOsr, 1)*[0:1:numViewSubCarriers-1]*fsub)).';
                figure(2);plot([0:1/viewSubCarrierOsr:numViewSubCarriers-1/viewSubCarrierOsr]*fsub, sincMat);xlabel('frequency (Hz)');ylabel('Amplitutde');title('Sub-carrier spacing and symbol rate relation');
                drawnow();
                firstTimePlotFlg = 1;
            end
            fsubSymRatioPrev = fsubSymRatio;
        else
            if viewSubCarriersPrev ~= viewSubCarriers
                close(figure(2));
            end
            firstTimePlotFlg = 0;
        end
        viewSubCarriersPrev = viewSubCarriers;
        
        
        % CHANNEL ESTIMATE
        if viewChEst3D == 1
            figure(3);
            surf(abs(chEstInterpArr));
            title('Channel Estimate - 3D view');zLimit = max(max(abs(chEstInterpArr)));xlabel('OFDM symbol');ylabel('sub-carrier');zlabel('Amplitude');
        else
            if viewChEst3D ~= viewChEst3DPrev
                close(figure(3));
            end
        end
        viewChEst3DPrev = viewChEst3D;
        
        if viewChEstMagPhase == 1
            figure(4);subplot(211);title('Channel Estimate for each OFDM symbol');
            plot(abs(chEstInterp), 'b-x');hold on;plot(abs(idealChanFreqDomain), 'r-o');hold off;xlabel('sub-carrier number');ylabel('Amplitude');grid on;legend('Estimated channel','Ideal channel');
            subplot(212);
            plot(angle(chEstInterp), 'b-x');hold on;plot(angle(idealChanFreqDomain), 'r-o');hold off;xlabel('sub-carrier number');ylabel('Phase');grid on;legend('Estimated channel','Ideal channel');
        else
            if viewChEstMagPhase ~= viewChEstMagPhasePrev
                close(figure(4));
            end
        end
        viewChEstMagPhasePrev = viewChEstMagPhase;
        
        if exitButton == 1;
            close(handles.figure1);
            return;
        end
    catch
        if exitButton == 1;
            close(handles.figure1);
            return;
        end
        set(handles.text18, 'String','Something has gone wrong! \n Exit and start over again');
        continue;
    end
    
    
end


function set_default_vals()

global bwSig;
global fs;
global fsChan;
global fsub;
global symRate;
global fdop;
global modTypeStr;
global pathGains;
global pathGainsTmp;
global pathDelays;
global pathDelaysTmp;
global Ncp;
global numOfdmSymsOnGrid;
global setDefaultFlag;
global viewTimeDomainWaveform;
global viewSubCarriers;
global pauseSim;
global exitButton;
global viewChEst3D;
global viewChEstMagPhase;
global useIdealChEstForEq;
global snrVal;
global cfoVal;
global phOffsetVal;
global viewEVM;

global bwSigNew;
global fsNew;
global fsChanNew;
global fsubNew;
global symRateNew;
global fdopNew;
global modTypeStrNew;
global pathGainsNew;
global pathDelaysNew;
global NcpNew;
global numOfdmSymsOnGridNew;
global snrValNew;
global cfoValNew;
global phOffsetValNew;

setDefaultFlag = 1;

bwSig = 1.4e6;
fs = 1.92e6;
fsChan = 0.01*fs;
fsub = 15e3;
symRate = fsub*1;
fdop = 300;
modTypeStr = 'QPSK';
pathGains = [0 -3 -6 -9];
pathGainsTmp = pathGains;
pathDelays = [0 0.4e-6 1e-6 1.5e-6]; % path delays in seconds
pathDelaysTmp = pathDelays;
Ncp = 20;
numOfdmSymsOnGrid = 14;
viewTimeDomainWaveform = 0;
viewSubCarriers = 0;
pauseSim = 0;
exitButton = 0;
viewChEst3D = 0;
viewChEstMagPhase = 0;
useIdealChEstForEq = 0;
snrVal = 1000;
cfoVal = 0;
phOffsetVal = 0;
viewEVM = 0;

bwSigNew	         = bwSig;
fsNew                = fs;
fsChanNew            = fsChan;
fsubNew              = fsub;
symRateNew           = symRate;
fdopNew              = fdop;
modTypeStrNew        = modTypeStr;
pathGainsNew         = pathGains;
pathDelaysNew        = pathDelays;
NcpNew               = Ncp;
numOfdmSymsOnGridNew = numOfdmSymsOnGrid;
snrValNew            = snrVal;
cfoValNew            = cfoVal;
phOffsetValNew       = phOffsetVal;



% EVM COMPUTATION
function [evmVal] = evm_compute(complexSyms, evmConstRef)
M = length(evmConstRef);
% find the reference point closest to the complex symbols
evmConstRefRepMat = (evmConstRef.'*ones(1, length(complexSyms)));
% find euclidean distance and compute EVM
[minDisVec, indexVec] = min((abs((ones(M, 1)*complexSyms) - evmConstRefRepMat).^2));
evmVal = sqrt(sum(minDisVec.^2)./sum(abs(evmConstRefRepMat([0:M:length(complexSyms)*M-1]+indexVec)).^2));
evmVal = round(evmVal*1e4)/1e4;

function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double
global pathDelaysNew;
try
    pathDelaysNew = eval(get(hObject,'String'));
catch
    set(handles.text18, 'String','enter valid path delays');
end

% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
global viewTimeDomainWaveform;
viewTimeDomainWaveform = get(hObject,'Value');

% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2

global viewSubCarriers
viewSubCarriers = get(hObject,'Value');


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pauseSim;
if pauseSim == 1
    pauseSim = 0;
else
    pauseSim = 1;
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global exitButton;
exitButton = 1;



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double
global cfoValNew;
try
    cfoValNew = eval(get(hObject,'String'));
catch
    set(handles.text18, 'String','enter a valid carrier frequency offset');
end


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
global viewChEst3D;
viewChEst3D = get(hObject,'Value');


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4
global useIdealChEstForEq;
useIdealChEstForEq = get(hObject,'Value');


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5

global viewChEstMagPhase;
viewChEstMagPhase = get(hObject,'Value');


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6
global viewEVM;
viewEVM = get(hObject, 'Value');


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global bwSig;
global fs;
global fsChan;
global fsub;
global symRate;
global fdop;
global modTypeStr;
global pathGains;
global pathDelays;
global Ncp;
global numOfdmSymsOnGrid;
global snrVal;
global cfoVal;
global phOffsetVal;

global bwSigNew;
global fsNew;
global fsChanNew;
global fsubNew;
global symRateNew;
global fdopNew;
global modTypeStrNew;
global pathGainsNew;
global pathDelaysNew;
global NcpNew;
global numOfdmSymsOnGridNew;
global snrValNew;
global cfoValNew;
global phOffsetValNew;

fsChan                 = fsChanNew;
modTypeStr             = modTypeStrNew;
Ncp                    = NcpNew;
numOfdmSymsOnGrid      = numOfdmSymsOnGridNew;
snrVal                 = snrValNew;
cfoVal                 = cfoValNew;
phOffsetVal 	        = phOffsetValNew;

set(handles.text18, 'String','');

if isequal(size(pathGainsNew), size(pathDelaysNew))
    if ~isequal(pathGainsNew, pathGains)
        set(handles.text18, 'String','path gains are updated');
    end
    pathGains = pathGainsNew;
    pathDelays = pathDelaysNew;
else
    set(handles.text18, 'String','the sizes of path gains and path delays are not compatible - consider revision');
    pathGainsStr = sprintf('%0.1f,', round(pathGains*10)/10);
    pathDelaysStr = sprintf('%0.1f,', round(pathDelays*1e6*10)/10);
    set(handles.edit7, 'String', ['[' pathGainsStr(1:end-1) ']']);
    set(handles.edit12, 'String', ['[' pathDelaysStr(1:end-1) ']*1e-6']);
    pathDelaysNew = pathDelays;
    pathGainsNew = pathGains;
end

if fsNew < bwSigNew
    set(handles.text18, 'String','choose sampling frequency greater than the bandwidth');
    set(handles.edit1, 'String', [num2str(fs/1e6) 'e6']);
    set(handles.edit2, 'String', [num2str(bwSig/1e6) 'e6']);
    bwSigNew = bwSig;
    fsNew = fs;
else
    bwSig                  = bwSigNew;
    fs                     = fsNew;
end

% Caution - DO NOT move the following condition from the current position
% bandwidth value needs to be frozen before sub-carrier spacing is determined
if fsubNew > (bwSig/4) % 4 is by choice
    if ~isequal(fsub, fsubNew)
        set(handles.text18, 'String','choose a smaller sub-carrier spacing');
        set(handles.edit3, 'String', [num2str(fsub/1e3) 'e3']);
        fsubNew = fsub;
    else
        set(handles.text18, 'String','bandwidth too low - forcing a higher bandwidth');
        set(handles.edit2, 'String', [num2str(fsub*4/1e3) 'e3']);
        bwSig = fsub*4;
        bwSigNew = bwSig;
    end
else
    fsub = fsubNew;
end

% Caution - DO NOT move the following condition from the current position
% bandwidth value needs to be frozen before symbol rate is determined
if symRateNew > (bwSig/4) % 4 is by choice
    set(handles.text18, 'String','choose a smaller OFDM symbol rate');
    set(handles.edit4, 'String', [num2str(symRate/1e3) 'e3']);
    symRateNew = symRate;
else
    symRate = symRateNew;
end

% Caution - DO NOT move the following condition from the current position
% Ncp has to be determined based on bandwidth and sub-carrier spacing - the
% order of the code matters
if NcpNew > (bwSig/fsub) % number of CP samples should be stricly less than the number of samples per OFDM symbol
    set(handles.text18, 'String','forcing a smaller cyclic prefix length');
    Ncp = floor(bwSig/fsub)-1;
    if Ncp < 0
        Ncp = 0;
    end
    set(handles.edit5, 'String', num2str(Ncp));drawnow();
    NcpNew = Ncp;
else
    Ncp = NcpNew;
end

if fdopNew < 5
    set(handles.text18, 'String','set doppler to at least 5 Hz');
    fdopNew = fdop;
else
    fdop = fdopNew;
end

if max(round(pathDelays*fs)) > round(bwSigNew/fsubNew)
    set(handles.text18, 'String','Path delays are too high. Forcing to only LOS path');
    pathDelays = 0;
    pathGains = 0;
    set(handles.edit7, 'String', [num2str(pathGains)]);
    set(handles.edit12, 'String', [num2str(pathDelays) '*1e-6']);
    pathDelaysNew = pathDelays;
    pathGainsNew = pathGains;
end
