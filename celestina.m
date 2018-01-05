function varargout = celestina(varargin)
% CELESTINA MATLAB code for celestina.fig
%      CELESTINA, by itself, creates a new CELESTINA or raises the existing
%      singleton*.
%
%      H = CELESTINA returns the handle to a new CELESTINA or the handle to
%      the existing singleton*.
%
%      CELESTINA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CELESTINA.M with the given input arguments.
%
%      CELESTINA('Property','Value',...) creates a new CELESTINA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before celestina_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to celestina_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help celestina

% Last Modified by GUIDE v2.5 23-Feb-2017 18:11:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @celestina_OpeningFcn, ...
                   'gui_OutputFcn',  @celestina_OutputFcn, ...
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


% --- Executes just before celestina is made visible.
function celestina_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to celestina (see VARARGIN)

% Choose default command line output for celestina
handles.output = hObject;
handles.ticks = [handles.waveforms_ax,handles.cross_corr, handles.time_ax,...
    handles.heat_a,handles.heat_b,handles.heat_c,handles.waveforms_ax,...
    handles.isi_c,handles.isi_a,handles.isi_b,handles.peak_distr,handles.cs_ix];
% Update handles structure
handles.black_axes = [handles.heat_a,handles.heat_b,handles.heat_c];

guidata(hObject, handles);
linkaxes([handles.heat_a,handles.heat_b,handles.heat_c],'xy');
linkaxes([handles.isi_a,handles.isi_b,handles.isi_c],'xy');
% UIWAIT makes celestina wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = celestina_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in load_gui_b.
function load_gui_b_Callback(hObject, eventdata, handles)
    h_figs = get(0,'children');
    h_fig = findobj(h_figs,'tag','wave_clus_figure');
    USER_DATA = get(h_fig,'userdata');
%     wc_handles = guidata(h_fig);
    if isempty(USER_DATA)
        set(handles.name_label,'String','Wave_clus GUI not found')
    end
    handles.classes = USER_DATA{6};
    handles.forced = USER_DATA{13};  
    handles.spikes = USER_DATA{2};
    handles.index = USER_DATA{3};
    handles.par = USER_DATA{1};
    set(handles.name_label,'String',[cd filesep USER_DATA{1}.nick_name])
    guidata(hObject,handles)
    update_figure(handles)

% hObject    handle to load_gui_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in ca_pm.
function ca_pm_Callback(hObject, eventdata, handles)
% hObject    handle to ca_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for i = handles.ticks
	cla(i)
end
legend(handles.waveforms_ax,'hide')
legend(handles.cs_ix,'hide')

% --- Executes during object creation, after setting all properties.
function ca_pm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ca_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


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


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in cb_pm.
function cb_pm_Callback(hObject, eventdata, handles)
% hObject    handle to cb_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for i = handles.ticks
	cla(i)
end
legend(handles.waveforms_ax,'hide')
legend(handles.cs_ix,'hide')

% --- Executes during object creation, after setting all properties.
function cb_pm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cb_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_b.
function plot_b_Callback(hObject, eventdata, handles)
    % hObject    handle to plot_b (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    contents = get(handles.ca_pm,'String'); 

    ca = str2num(contents{get(handles.ca_pm,'Value')});
    cb = str2num(contents{get(handles.cb_pm,'Value')});
    handles.ca = ca;
    handles.cb = cb;
    
    sp_a = handles.classes == ca;
    sp_b = handles.classes == cb;
    sp_c = handles.classes  == ca | handles.classes == cb;
    sps = handles.spikes;
    %waveform axes
    cla(handles.waveforms_ax)
    hold(handles.waveforms_ax,'on')
    av = mean(sps(sp_a,:),1);
    sd = std(sps(sp_a,:),1,1);
    p(1) = plot(handles.waveforms_ax,av,'-r','linewidth',1);
    plot(handles.waveforms_ax,av+sd,':r','linewidth',1);
    plot(handles.waveforms_ax,av-sd,':r','linewidth',1);
    av = mean(sps(sp_b,:),1);
    sd = std(sps(sp_b,:),1,1);
    p(2) = plot(handles.waveforms_ax,av,'-b','linewidth',1);
    plot(handles.waveforms_ax,av+sd,':b','linewidth',1);
    plot(handles.waveforms_ax,av-sd,':b','linewidth',1);
    av = mean(sps(sp_c,:),1);
    sd = std(sps(sp_c,:),1,1);
    p(3) = plot(handles.waveforms_ax,av,'-k','linewidth',1);
    plot(handles.waveforms_ax,av+sd,':k','linewidth',1);
    plot(handles.waveforms_ax,av-sd,':k','linewidth',1);
    legend(handles.waveforms_ax,p,['C:' num2str(ca)],['C:' num2str(cb)],'Merge','Location','Best')
    xlim(handles.waveforms_ax,[1 length(av)]);

    %heat maps
    plot_heatmaps(handles,sps,ca,cb,sp_a,sp_b,sp_c);

    %isis
    nbins= 100;
    bin_step = 1;
    ix = handles.index; %to ms
    % Calculates # ISIs < 3ms  
    times = diff(ix(sp_a)); 
    ra = nnz(times < 3); 
    [N,X]=hist(times,0:bin_step:nbins);
    bar(handles.isi_a,X(1:end-1),N(1:end-1),'r')
    
    times = diff(ix(sp_b)); 
    rb = nnz(times < 3); 
    [N,X]=hist(times,0:bin_step:nbins);
    bar(handles.isi_b,X(1:end-1),N(1:end-1),'b')
   
    times = diff(ix(sp_c)); 
    rc = nnz(times < 3); 
    [N,X]=hist(times,0:bin_step:nbins);
    bar(handles.isi_c,X(1:end-1),N(1:end-1),'k')
    nmax = max(N(1:end-1));
    if nmax>0
        ylim(handles.isi_c,[0 nmax]);
    end
    xlim(handles.isi_c,[0 nbins]);
    xlabel(handles.isi_c,'ISI (ms)');

    text(nbins,nmax,[num2str(rc) ' in < 3ms '],'Parent',handles.isi_c,'HorizontalAlignment','right','VerticalAlignment','top')
    text(nbins,nmax,[num2str(rb) ' in < 3ms '],'Parent',handles.isi_b,'HorizontalAlignment','right','VerticalAlignment','top')
    text(nbins,nmax,[num2str(ra) ' in < 3ms '],'Parent',handles.isi_a,'HorizontalAlignment','right','VerticalAlignment','top')
    %time plot
    plot_time(handles,sps,sp_a,sp_b,sp_c,ix)
        
    %cross correlation
    Na = hist(ix(sp_a),min(ix(sp_c)):10:max(ix(sp_c)));
    Nb = hist(ix(sp_b),min(ix(sp_c)):10:max(ix(sp_c)));
    maxlag = 400;
    [acor,lag] = xcorr(Na,Nb,maxlag);
    bar(handles.cross_corr,lag,acor,'FaceColor',[0.5,0.25,0.6])
    xlim(handles.cross_corr,[-maxlag maxlag])
    xlabel(handles.cross_corr,'time(ms)')
    ylabel(handles.cross_corr,'Cross-Correlation' ,'fontunits',   'normalized', 'FontSize', 0.08)
    yl = ylim(handles.cross_corr);
    ylim(handles.cross_corr,[0 yl(2)]);
    
    %peak distribution
    if strcmp(handles.par.detection,'both') || any(isnan(handles.par.detection))
    	peak = @(x) max(abs(x),[],2);
    elseif strcmp(handles.par.detection,'neg')
        peak = @(x) max(x,[],2);
    elseif strcmp(handles.par.detection,'pos')
        peak = @(x) min(x,[],2);
    end
    apeak = peak(sps(sp_a,:));
    bpeak = peak(sps(sp_b,:));
    max_ap = max(apeak);
    min_ap = min(apeak);
    max_bp = max(bpeak);
    min_bp = min(bpeak);
    hd = (max(max_ap,max_bp) - min(min_ap,min_bp))/100;
    [Pa xa] = hist(apeak,min_ap:hd:max_ap);
    [Pb xb]= hist(bpeak,min_bp:hd:max_bp);
    hold(handles.peak_distr,'on')
    bar(handles.peak_distr,xa,Pa,'r')
    bar(handles.peak_distr,xb,Pb,'b')
    ylabel(handles.peak_distr,'Distrib. of amplitudes' ,'fontunits',   'normalized', 'FontSize', 0.08)
%     set(handles.peak_distr, 'Color', 'none')
    xlabel(handles.peak_distr,'amplitude')
    xlim(handles.peak_distr,'auto')
    %cumulative sum
    cla(handles.cs_ix)
    hold(handles.cs_ix,'on')
    plot(handles.cs_ix,ix(sp_a),1:nnz(sp_a),'r')
    plot(handles.cs_ix,ix(sp_b),1:nnz(sp_b),'b')
    plot(handles.cs_ix,ix(sp_c),1:nnz(sp_c),'k')
    xlim(handles.cs_ix,[0 max(ix(sp_c))])
    ylim(handles.cs_ix,'auto')
    ylim(handles.cs_ix,ylim(handles.cs_ix).*[0 1])
    legend(handles.cs_ix,['C:' num2str(ca)],['C:' num2str(cb)],'Merge','Location','Best')
    ylabel(handles.cs_ix,'cumulat. spike count','fontunits',   'normalized', 'FontSize', 0.08)
    xlabel(handles.cs_ix,'time(ms)')
%     set(handles.cs_ix, 'Color', 'none')
    
    for i = handles.ticks
        set(i,  'fontunits',   'normalized', 'FontSize', 0.06)
    end
    guidata(hObject,handles)
    
function heat_scale_pm_Callback(hObject, eventdata, handles)
    if isempty(handles.ca)
        return
    end
    sp_a = handles.classes == handles.ca;
    sp_b = handles.classes == handles.cb;
    sp_c = handles.classes  == handles.ca | handles.classes == handles.cb;
    sps = handles.spikes;
    plot_heatmaps(handles,sps,handles.ca,handles.cb,sp_a,sp_b,sp_c)

% --- Executes during object creation, after setting all properties.
function heat_scale_pm_CreateFcn(hObject, eventdata, handles)
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in time_plot_pm.
function time_plot_pm_Callback(hObject, eventdata, handles)
    if isempty(handles.ca)
        return
    end
    sp_a = handles.classes == handles.ca;
    sp_b = handles.classes == handles.cb;
    sp_c = handles.classes  == handles.ca | handles.classes == handles.cb;
    plot_time(handles,handles.spikes,sp_a,sp_b,sp_c,handles.index)

    

% --- Executes during object creation, after setting all properties.
function time_plot_pm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_plot_pm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function update_figure(h)
    cl = unique(h.classes);
    cl = cl(cl>0);
    for i=1:length(cl)
        text{i} = num2str(cl(i));
    end
    set(h.ca_pm,'String',text);
    set(h.cb_pm,'String',text);

    for i = h.ticks
        cla(i)
    end
    legend(h.waveforms_ax,'hide')
    legend(h.cs_ix,'hide')
    handles.ca = [];
    handles.cb = [];
% --- Executes on button press in load_b.
function load_b_Callback(hObject, eventdata, handles)
% hObject    handle to load_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [filename, pathname] = uigetfile('times_*.mat','Select file');
    if isempty(filename)
        return
    end
    load([pathname filesep filename],'cluster_class','spikes','forced','par')
    handles.classes = cluster_class(:,1);
    if exist('force','var')
        handles.forced = forced;  
    else
        handles.forced = [];  
    end
    handles.spikes = spikes;
    handles.index = cluster_class(:,2);
    handles.par = par;
    set(handles.name_label,'String',[pathname filesep filename])
    guidata(hObject,handles)
    update_figure(handles)



function isi_step_et_Callback(hObject, eventdata, handles)
% hObject    handle to isi_step_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of isi_step_et as text
%        str2double(get(hObject,'String')) returns contents of isi_step_et as a double


% --- Executes during object creation, after setting all properties.
function isi_step_et_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isi_step_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function isi_max_et_Callback(hObject, eventdata, handles)
% hObject    handle to isi_max_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of isi_max_et as text
%        str2double(get(hObject,'String')) returns contents of isi_max_et as a double


% --- Executes during object creation, after setting all properties.
function isi_max_et_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isi_max_et (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plot_time(handles,sps,sp_a,sp_b,sp_c,ix)
    contents = cellstr(get(handles.time_plot_pm,'String'));% returns heat_scale_pm contents as cell array
    valuetype = contents{get(handles.time_plot_pm,'Value')}; %returns selected item from heat_scale_pm
    
    cla(handles.time_ax)
    hold(handles.time_ax,'on')
    
   switch valuetype
   case 'min'
        apeak = min(sps(sp_a,:),[],2);
        bpeak = min(sps(sp_b,:),[],2);
   case 'max'
        apeak = max(sps(sp_a,:),[],2);
        bpeak = max(sps(sp_b,:),[],2);
   case 'ptp'
        apeak = max(sps(sp_a,:),[],2)-min(sps(sp_a,:),[],2);
        bpeak = max(sps(sp_b,:),[],2)-min(sps(sp_b,:),[],2);       
   end
    
    plot(handles.time_ax,ix(sp_a),apeak,'.r', 'markersize', 5)
    plot(handles.time_ax,ix(sp_b),bpeak,'.b', 'markersize', 5)
    xlabel(handles.time_ax,'time(ms)','fontunits',   'normalized', 'FontSize', 0.08);
    xlim(handles.time_ax,[min(ix(sp_c)), max(ix(sp_c))]);
    
function plot_heatmaps(handles,sps,ca,cb,sp_a,sp_b,sp_c)
%heat maps
    contents = cellstr(get(handles.heat_scale_pm,'String'));% returns heat_scale_pm contents as cell array
    islog = strcmp(contents{get(handles.heat_scale_pm,'Value')},'Log'); %returns selected item from heat_scale_pm
    
    lsp = size(sps,2);
    x = repmat((1:lsp)',nnz(sp_a),1);
    y = reshape(sps(sp_a,:)',[nnz(sp_a)*lsp,1]);
    [aux,aux_c] = hist3([x,y],[lsp lsp]);
    if islog
        aux=log10(aux+1);
    end
    colormap(handles.heat_a,'hot')
    pcolor(handles.heat_a,aux_c{1},aux_c{2},aux')
    shading(handles.heat_a,'Flat');%interp
    
    x = repmat((1:lsp)',nnz(sp_b),1);
    y = reshape(sps(sp_b,:)',[nnz(sp_b)*lsp,1]);
    [aux,aux_c] = hist3([x,y],[lsp lsp]);
    if islog
        aux=log10(aux+1);
    end
    colormap(handles.heat_b,'hot')
    pcolor(handles.heat_b,aux_c{1},aux_c{2},aux')
    shading(handles.heat_b,'Flat');%interp
    x = repmat((1:lsp)',nnz(sp_c),1);
    y = reshape(sps(sp_c,:)',[nnz(sp_c)*lsp,1]);
    [aux,aux_c] = hist3([x,y],[lsp lsp]);
    if islog
        aux=log10(aux+1);
    end
    colormap(handles.heat_c,'hot')
    pcolor(handles.heat_c,aux_c{1},aux_c{2},aux')
    shading(handles.heat_c,'Flat');%interp
    xlabel(handles.heat_c,'samples')
    
    ylabel(handles.heat_a, ['C:' num2str(ca) ' (#' num2str(nnz(sp_a)) ')'],'fontunits',   'normalized', 'FontSize', 0.10)
    ylabel(handles.heat_b, ['C:' num2str(cb) ' (#' num2str(nnz(sp_b)) ')'],  'fontunits',   'normalized', 'FontSize', 0.10)
    ylabel(handles.heat_c, ['Merge (#' num2str(nnz(sp_c)) ')'], 'FontSize', 12)
    set(handles.black_axes,'Color',[0 0 0 ]);

