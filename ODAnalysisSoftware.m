function varargout = ODAnalysisSoftware(varargin)
% ODANALYSISSOFTWARE MATLAB code for ODAnalysisSoftware.fig
%      ODANALYSISSOFTWARE, by itself, creates a new ODANALYSISSOFTWARE or raises the existing
%      singleton*.
%
%      H = ODANALYSISSOFTWARE returns the handle to a new ODANALYSISSOFTWARE or the handle to
%      the existing singleton*.
%
%      ODANALYSISSOFTWARE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ODANALYSISSOFTWARE.M with the given input arguments.
%
%      ODANALYSISSOFTWARE('Property','Value',...) creates a new ODANALYSISSOFTWARE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ODAnalysisSoftware_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ODAnalysisSoftware_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ODAnalysisSoftware

% Last Modified by GUIDE v2.5 13-Aug-2019 13:49:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ODAnalysisSoftware_OpeningFcn, ...
    'gui_OutputFcn',  @ODAnalysisSoftware_OutputFcn, ...
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


% --- Executes just before ODAnalysisSoftware is made visible.
function ODAnalysisSoftware_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ODAnalysisSoftware (see VARARGIN)

% Choose default command line output for ODAnalysisSoftware
currentFolder = pwd ; 
addpath([currentFolder '\Data'])
addpath([currentFolder '\functions'])
addpath([currentFolder '\Image'])


baseFileName = 'sunyopt.jpg';
rgbImage = imread(baseFileName);
axes(handles.axesLogo);
imshow(rgbImage, []);




handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ODAnalysisSoftware wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ODAnalysisSoftware_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbuttonImageSelection.
function pushbuttonImageSelection_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonImageSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Condition = [get(handles.radiobuttonHuman,'value') get(handles.radiobuttonMacaque,'value')...
    get(handles.radiobuttonCat,'value') get(handles.radiobuttonOther,'value')];

if Condition(1) == 1
    
    load('workspaceHumanfig3.mat')
    surrounding_exclude_scalebar = surrounding;
    %surrounding_exclude_scalebar(1056:1062,65:281) = 0;
    od_human_for_video = 2*input_bw.*v1_region + surrounding_exclude_scalebar*1;
    
    
    cla(handles.axesImage,'reset')
    %cla(handles.axesResult,'reset'),axis off
    %cla(handles.axesOrientation,'reset'),axis off
    axes(handles.axesImage);
    imagesc(od_human_for_video)
    %colormap([0 0 0;1 1 1;0.7 0.7 0.7])
    colormap([0 0 0;0.7 0.7 0.7;1 1 1])
    axis off
    %text(60, 1010,'10 mm','fontsize',20)
    title('Ocular Dominance of Human','FontSize',16)
    
    set(handles.pushbuttonFeatureSelction,'visible','on')
    set(handles.pushbuttonASF,'visible','on')
    set(handles.pushbuttonMakeVideo,'visible','on')
    
elseif Condition(2) == 1
    
    load('workspaceMacaque.mat')
    ODImage = input_bw.*v1_region*2 + surrounding;
    cla(handles.axesImage,'reset')
    axes(handles.axesImage)
    imagesc(ODImage)
    colormap([0 0 0;0.7 0.7 0.7;1 1 1])
    axis off
    title('Ocular Dominance of Macaque','FontSize',16)
    
    set(handles.pushbuttonFeatureSelction,'visible','on')
    set(handles.pushbuttonASF,'visible','on')
    set(handles.pushbuttonMakeVideo,'visible','on')
    
elseif Condition(3) == 1
    
    load('workspaceCat.mat')
    ODImage = input_bw.*v1_region*2 + surrounding;
    cla(handles.axesImage,'reset')
    axes(handles.axesImage)
    imagesc(ODImage)
    colormap([0 0 0;0.7 0.7 0.7;1 1 1])
    axis off
    title('Ocular Dominance of Cat','FontSize',16)
    
    set(handles.pushbuttonFeatureSelction,'visible','on')
    set(handles.pushbuttonASF,'visible','on')
    set(handles.pushbuttonMakeVideo,'visible','on')
    
elseif Condition(4) == 1
    
    set(handles.pushbuttonFeatureSelction,'visible','on')
    set(handles.pushbuttonASF,'visible','off')
    set(handles.pushbuttonMakeVideo,'visible','off')
    
    [p,q] = uigetfile('*.*','Select the data file');
    InputImage = imread([q,p]);
    handles.InputImage = InputImage;
    
    cla(handles.axesImage,'reset')
    axes(handles.axesImage)
    imshow(InputImage)
    
end


handles.Condition = Condition;
guidata(hObject,handles);


% --- Executes on button press in pushbuttonFeatureSelction.
function pushbuttonFeatureSelction_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonFeatureSelction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Condition = handles.Condition;


axes(handles.axesImage);
[col,row] = ginput(1);
col = round(col);
row = round(row);

handles.SelectedPoint = [row col];

if Condition(1) == 1
    
    load('workspaceHumanfig3.mat')
    surrounding_exclude_scalebar = surrounding;
    %surrounding_exclude_scalebar(1056:1062,65:281) = 0;
    ODImage = input_bw.*v1_region*2 + surrounding_exclude_scalebar;
    
elseif Condition(2) == 1
    
    load('workspaceMacaque.mat')
    ODImage = input_bw.*v1_region*2 + surrounding;
    
elseif Condition(3) == 1
    
    load('workspaceCat.mat')
    ODImage = input_bw.*v1_region*2 + surrounding;
    
elseif Condition(4) == 1
    
    %functionOtherImageAnalysis(handles)
    
end

%   figure,imshow(surrounding)
%   figure,imshow(ODImage)

if Condition(4) ~= 1
    step_hist  = 30;
    initial_edge = -(step_hist/2);
    last_edge = 360 + initial_edge;
    edge_range = (initial_edge:step_hist:last_edge);
    bin_plot= (0:step_hist:360) *(pi/180);
    
    
    [L_ipsi,N_region_ipsi] = bwlabel(~input_bw.*v1_region);%ipsi black
    [L_contra,N_region_contra] = bwlabel(input_bw.*v1_region);%contra _ white
    %   figure,imagesc(L_contra)
    %   figure,imagesc(v1_region)
    
    if L_ipsi(row,col) > 0
        i_ipsi  = L_ipsi(row,col);
        selected_region = L_ipsi  == i_ipsi ;
        angle_line = output_orientation_ipsi{i_ipsi};
        thick = output_thickness_ipsi{i_ipsi};
        
        cla(handles.axesImage,'reset')
        axes(handles.axesImage); 
        imagesc(ODImage)
        colormap([0 0 0;0.7 0.7 0.7;1 1 1])
        axis off
        %text(60, 1010,'10 mm','fontsize',20)
        title(sprintf(' Stripe Length = %.0f pixels (%.2f mm)',npoint_ipsi(i_ipsi),npoint_ipsi(i_ipsi)*pixel2um/1000),'fontsize',18)
        [selected_region_row,selected_region_col] = find(edge(selected_region));
        hold on, plot(selected_region_col,selected_region_row,'.r')
        
        cla(handles.axesResult)
        set(handles.axesResult,'visible','on')
        axes(handles.axesResult)
        thick(thick<0) = 0;
        [N_thick_ipsi,plot2] = hist_mid_line(thick,nbins_thickness_ipsi,'k');
        axis([.9*min(thick) 1.1*max(thick) -.05 1])
        xlabel('Stripe Width ','fontsize',14)
        ylabel('Frequency','fontsize',14)
        title(sprintf('                                          Stripe Width  %.2f pixels (%.3f mm)',thickness_ipsi(i_ipsi),thickness_ipsi(i_ipsi)*pixel2um/1000),'fontsize',16)
        ax_width = gca;
        get( ax_width );
        set( ax_width, 'Color', [0.7,0.7,0.7] )
        ax_width.TickDir = 'out';
        ax_width.TickLength = [0.02 0.02];
        box(ax_width,'off')
        
        axes(handles.axesOrientation)
        set(handles.axesOrientation,'visible','on')
        angle_line_0_180 = angle_line .* (angle_line>0) + (angle_line+180) .* (angle_line<0);
        angle_line_full_range_I = cat(2,angle_line_0_180,angle_line_0_180+180);
        [frequency_ipsi,bin_edge_ipsi] = histcounts(angle_line_full_range_I * (pi/180) ,edge_range * (pi/180),'Normalization','Probability'); %edge is the bar start and end points
        frequency_ipsi(1) = frequency_ipsi(floor(end/2)+1);% 0 and 180 should be equal
        frequency_ipsi(end+1) = frequency_ipsi(1); %    for smooth plotting
        polar_histogram_predefined_edge_makevideo(bin_plot,frequency_ipsi,max(frequency_ipsi),'-',[0 0 0]/255);
        xlabel('Stripe Angle','fontsize',14)
        
    elseif L_contra(row,col) > 0
        i_contra  = L_contra(row,col);
        selected_region = L_contra  == i_contra ;
        angle_line = output_orientation_contra{i_contra};
        thick = output_thickness_contra{i_contra};
        
        cla(handles.axesImage,'reset')
        axes(handles.axesImage)
        imagesc(ODImage)
        colormap([0 0 0;0.7 0.7 0.7;1 1 1])
        axis off
        %text(60, 1010,'10 mm','fontsize',20)
        title(sprintf('Stripe Length = %.0f pixels (%.2f mm)',npoint_contra(i_contra),npoint_contra(i_contra)*pixel2um/1000),'fontsize',18)
        [selected_region_row,selected_region_col] = find(edge(selected_region));
        hold on, plot(selected_region_col,selected_region_row,'.r')
        
        cla(handles.axesResult)
        axes(handles.axesResult)
        set(handles.axesResult,'visible','on')
        thick(thick<0) = 0;
        [N_thick_contra,plot2] = hist_mid_line(thick,nbins_thickness_contra,'w');
        axis([.9*min(thick) 1.1*max(thick) -.05 1])
        xlabel('Stripe Width ','fontsize',14)
        ylabel('Frequency','fontsize',14)
        title(sprintf('                                           Stripe Width = %.2f pixels (%.3f mm)',thickness_contra(i_contra),thickness_contra(i_contra)*pixel2um/1000),'fontsize',16)
        ax_width = gca;
        get( ax_width );
        set( ax_width, 'Color', [0.7,0.7,0.7] )
        ax_width.TickDir = 'out';
        ax_width.TickLength = [0.02 0.02];
        box(ax_width,'off')
        
        axes(handles.axesOrientation)
        set(handles.axesOrientation,'visible','on')
        angle_line_0_180 = angle_line .* (angle_line>0) + (angle_line+180) .* (angle_line<0);
        angle_line_full_range_C = cat(2,angle_line_0_180,angle_line_0_180+180);
        [frequency_contra,bin_edge_contra] = histcounts(angle_line_full_range_C * (pi/180) ,edge_range * (pi/180),'Normalization','Probability'); %edge is the bar start and end points
        frequency_contra(1) = frequency_contra(floor(end/2)+1);% 0 and 180 should be equal
        frequency_contra(end+1) = frequency_contra(1); %    for smooth plotting
        polar_histogram_predefined_edge_makevideo(bin_plot,frequency_contra,max(frequency_contra),'-',[1 1 1]);
        xlabel('Stripe Angle','fontsize',14)
    end
    
end
guidata(hObject,handles);


% --- Executes on button press in pushbuttonASF.
function pushbuttonASF_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonASF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Condition = handles.Condition;
cla(handles.axesResult,'reset')
cla(handles.axesOrientation,'reset')
set(handles.axesResult,'visible','off')
set(handles.axesOrientation,'visible','off')
functionDipoleFourierSoftware(Condition,handles)
guidata(hObject,handles);



% --- Executes on button press in pushbuttonMakeVideo.
function pushbuttonMakeVideo_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonMakeVideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

functionMakeMovieSoftware(handles)

guidata(hObject,handles);
