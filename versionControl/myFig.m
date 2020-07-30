function varargout = myFig(varargin)
% MYFIG MATLAB code for myFig.fig
%      MYFIG, by itself, creates a new MYFIG or raises the existing
%      singleton*.
%
%      H = MYFIG returns the handle to a new MYFIG or the handle to
%      the existing singleton*.
%
%      MYFIG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MYFIG.M with the given input arguments.
%
%      MYFIG('Property','Value',...) creates a new MYFIG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before myFig_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to myFig_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help myFig

% Last Modified by GUIDE v2.5 04-May-2018 11:36:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @myFig_OpeningFcn, ...
                   'gui_OutputFcn',  @myFig_OutputFcn, ...
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


% --- Executes just before myFig is made visible.
function myFig_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to myFig (see VARARGIN)

% Choose default command line output for myFig
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes myFig wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = myFig_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in playAudio.
function playAudio_Callback(hObject, eventdata, handles)
% hObject    handle to playAudio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global x
global fs
% global rfinal
soundsc(x,fs);

% --------------------------------------------------------------------
function openFile_Callback(hObject, eventdata, handles)
% hObject    handle to openFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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

% plot original signal
set(handles.axes1,'HandleVisibility','ON');
axes(handles.axes1);
axis on;
% plot(t3,x/(max(x)),'Color',[0.3137 0.3137 0.3137]);
plot(t3,x,'Color',[0.3137 0.3137 0.3137]);
ylim([-2 2]);
grid on;
clc;


% --- Executes on button press in Spectrogram.
function Spectrogram_Callback(hObject, eventdata, handles)
% hObject    handle to Spectrogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% calculate spectrogram
global x
global fs

[magSpec,~,~] = spectrogram(x,hamming(512),256,fs);

dbSpec = mag2db(abs(magSpec));

% plot original signal
set(handles.axes2,'HandleVisibility','ON');
axes(handles.axes2);
axis on;
imagesc(dbSpec); axis xy;
% ylim([-2 2]);
grid on;
clc;


% --- Executes on button press in featureExtraction.
function featureExtraction_Callback(hObject, eventdata, handles)
% hObject    handle to featureExtraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global x
global fs

% calcualte MFCCs

% calculate LFCCs


% calculate temporal and frequency features


% --- Executes on button press in showFeatures.
function showFeatures_Callback(hObject, eventdata, handles)
% hObject    handle to showFeatures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
