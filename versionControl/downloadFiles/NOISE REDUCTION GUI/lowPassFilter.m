function varargout = lowPassFilter(varargin)
 
% Begin initialization code
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lowPassFilter_OpeningFcn, ...
                   'gui_OutputFcn',  @lowPassFilter_OutputFcn, ...
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

% End initialization code  
% Executes just before lowPassFilter is made visible.
% --------------------------------------------------------------------
 
function file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data set(handles.axes8,'HandleVisibility','OFF');

% get file from user
clear global;
global FileName
global PathName
[FileName,PathName] = uigetfile('*.wav','Select the Wav-file');

% reading wav file
global x
global fs
[x,fs]=audioread(fullfile(PathName, FileName));
global t3
t3=(0:length(x)-1)/fs;

% plotting original signal
set(handles.axes16,'HandleVisibility','ON');
axes(handles.axes16);
axis on;
plot(t3,x/(max(x)),'Color',[0.3137 0.3137 0.3137]);
grid on;
clc;

function lowPassFilter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data 
% varargin   command line arguments to lowPassFilter  
% Choose default command line output for lowPassFilter
handles.output = hObject;
 
% Update handles structure
guidata(hObject, handles);
 
% UIWAIT makes lowPassFilter wait for user response (see UIRESUME)
% uiwait(handles.figure1);
 
function order_Callback(hObject, eventdata, handles)
% hObject    handle to order 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data  
% --- Executes during object creation, after setting all properties.

function order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
% hObject    handle to ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global u
global x
global fs
t=(0:length(x)-1)/fs;
L=length(x);

% fourier transform of signal
NFFT = 2^nextpow2(L);
Y=fft(u,NFFT);

% plotting frequecy spectrum
f = fs/2*linspace(0,1,NFFT/2);
om22 = str2double(get(handles.order,'String'));
be = str2double(get(handles.db,'string'));
fq = str2double(get(handles.cutoff,'string'));
be2=10^(be/20);
Fs2=length(f);
temp1=(fq*Fs2)/max(f);
temp2=(pi*temp1)/Fs2;
yelpha=(cos(pi/(2*om22)))/((cos(temp2/2))*(cosh((1/om22)*acosh(be2))));
set(handles.alpha,'String',yelpha)
fnew=20*log10(abs(Y(1:length(f))));
set(handles.axes6,'HandleVisibility','ON');
axes(handles.axes6);
plot(f,fnew,'Color',[0.3137 0.3137 0.3137]);
grid on;

%Filter Design Starts Here
b=be2;
m=om22;
a=[1:m];
i=1;
freq=0;
Fs2=length(f);
 
h=1;
p=1;
omegak=((2*a-1)*pi)/(2*m);
omegam=2*acos((cos(omegak))/(((yelpha)*(cosh(1/m*acosh(b))))));
while(freq<Fs2)
   om=2*pi*freq/(2*Fs2);
   k=exp(j*om);
   while(i<=m)
      h=h*(k-exp(j*omegam(i)));
      i=i+1;
   end
   h1(p)=abs(h);
   p=p+1;
   freq=freq+1;
   h=1;
   i=1;
end
h2=(h1/max(h1));
h3=20*log10(h2);
set(handles.axes11,'HandleVisibility','ON');
axes(handles.axes11);
axis tight;
plot(f,h3,'Color',[0.3137 0.3137 0.3137]);
grid on;
mul=(Y(1:length(f))).*h2';
mul2=abs(mul);
mulnew=20*log10(mul2);
set(handles.axes13,'HandleVisibility','ON');
axes(handles.axes13);
axis tight;
plot(f,mulnew,'Color',[0.3137 0.3137 0.3137]);
grid on;
mul3=ifft(mul,NFFT);
mul4=abs(mul3);
mul5=max(mul4);
mul6=mul4/mul5;
 
rmul=real(mul3(1:length(x)));
set(handles.axes12,'HandleVisibility','ON');
axes(handles.axes12);
axis tight;
plot(t,rmul/(max(rmul)),'Color',[0.3137 0.3137 0.3137]);
grid on;
global rfinal
rfinal=rmul/(max(rmul));

% storing the final filtered signal
audiowrite('ksfinal.wav', rfinal(1:length(x)), fs);
clc;


% --- Outputs from this function are returned to the command line.
function varargout = lowPassFilter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data  
% Get default command line output from handles structure
varargout{1} = handles.output;
 
function db_Callback(hObject, eventdata, handles)
% hObject    handle to db 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data
 
 
% --- Executes during object creation, after setting all properties.
function db_CreateFcn(hObject, eventdata, handles)
% hObject    handle to db 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
function cutoff_Callback(hObject, eventdata, handles)
% hObject    handle to cutoff
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data 

% --- Executes during object creation, after setting all properties.
function cutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cutoff 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
function alpha_Callback(hObject, eventdata, handles)
% hObject    handle to alpha 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data 

 
% --- Executes during object creation, after setting all properties.
function alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
function inp_Callback(hObject, eventdata, handles)
% hObject    handle to inp 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data  

% --- Executes during object creation, after setting all properties.
function inp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inp 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 
% --- Executes on button press in play1.
function play1_Callback(hObject, eventdata, handles)
% hObject    handle to play1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data 
global x
global fs
global rfinal
audioplayer(rfinal(1:length(x)),fs);
 
% --- Executes on button press in play2.
function play2_Callback(hObject, eventdata, handles)
% hObject    handle to play2
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data 
global u
global fs
audioplayer(u,fs);
 
% --- Executes on button press in original.
function original_Callback(hObject, eventdata, handles)
% hObject    handle to original 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data 
global x
global fs
audioplayer(x,fs);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global x
white = 0.1*rand(1,(length(x)));
NFFT2 = 2^nextpow2(length(white));
Y2=fft(white,NFFT2);
global fs
f2 = fs/2*linspace(0,1,NFFT2/2);
fnew2=20*log10(abs(Y2(1:length(f2))));
fnew2(1:(length(fnew2)/2))=0;
q1=ifft(Y2,NFFT2);
q2=abs(q1);
q3=max(q2);
q4=q2/q3;
global q5
q5=0.1*(real(q4(1:length(x))));
global fs
pin=str2double(get(handles.freqn,'String'));
cpnoise =0.1* cos(2*pi*pin*(0:length(x)-1)/fs)';
global u
u = x;
if (get(handles.pushbutton5,'Value') == get(hObject,'Max'))
    u = u + q5';
end
if (get(handles.pushbutton6,'Value') == get(hObject,'Max'))
    u = u + cpnoise;
end

global fs

% saving the noisy file
audiowrite('mynoisyfile.wav', u, fs);
t=(0:length(u)-1)/fs;
set(handles.axes8,'HandleVisibility','ON');
axes(handles.axes8);
axis on;
plot(t,u/(max(u)),'Color',[0.3137 0.3137 0.3137]);
grid on;
% plotting frequecy spectrum
t=(0:length(x)-1)/fs;
L=length(x);

% fourier transform of signal
NFFT = 2^nextpow2(L);
Y=fft(u,NFFT);
f = fs/2*linspace(0,1,NFFT/2);
om22 = str2double(get(handles.order,'String'));
be = str2double(get(handles.db,'string'));
fq = str2double(get(handles.cutoff,'string'));
be2=10^(be/20);
Fs2=length(f);
temp1=(fq*Fs2)/max(f);
temp2=(pi*temp1)/Fs2;
yelpha=(cos(pi/(2*om22)))/((cos(temp2/2))*(cosh((1/om22)*acosh(be2))));
set(handles.alpha,'String',yelpha)
fnew=20*log10(abs(Y(1:length(f))));

set(handles.axes6,'HandleVisibility','ON');
axes(handles.axes6);
plot(f,fnew,'Color',[0.3137 0.3137 0.3137]);
grid on;

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global x
global fs
pin=str2double(get(handles.freqn,'String'));
cpnoise =0.1* cos(2*pi*pin*(0:length(x)-1)/fs)';
global u
u = x;
length(u)
length(cpnoise)
if (get(handles.pushbutton6,'Value') == get(hObject,'Max'))
    u = x + cpnoise;
end
if (get(handles.pushbutton5,'Value') == get(hObject,'Max'))
    global q5
    u = u + q5';
end
global fs

% plotting frequecy spectrum
t=(0:length(x)-1)/fs;
L=length(x);

% fourier transform of signal
NFFT = 2^nextpow2(L);
Y=fft(u,NFFT);
f = fs/2*linspace(0,1,NFFT/2);
om22 = str2double(get(handles.order,'String'));
be = str2double(get(handles.db,'string'));
fq = str2double(get(handles.cutoff,'string'));
be2=10^(be/20);
Fs2=length(f);
temp1=(fq*Fs2)/max(f);
temp2=(pi*temp1)/Fs2;
yelpha=(cos(pi/(2*om22)))/((cos(temp2/2))*(cosh((1/om22)*acosh(be2))));
set(handles.alpha,'String',yelpha)
fnew=20*log10(abs(Y(1:length(f))));
set(handles.axes6,'HandleVisibility','ON');
axes(handles.axes6);
plot(f,fnew,'Color',[0.3137 0.3137 0.3137]);
grid on;

% saving the noisy file
audiowrite('mynoisyfile.wav',u,fs);
t=(0:length(u)-1)/fs;
set(handles.axes8,'HandleVisibility','ON');
axes(handles.axes8);
axis on;
plot(t,u/(max(u)),'Color',[0.3137 0.3137 0.3137]);
grid on;
 


function freqn_Callback(hObject, eventdata, handles)
% hObject    handle to freqn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freqn as text
%        str2double(get(hObject,'String')) returns contents of freqn as a double


% --- Executes during object creation, after setting all properties.
function freqn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freqn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


