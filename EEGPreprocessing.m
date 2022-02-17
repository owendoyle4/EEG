%% Preprocessing steps
% 1. initialize: enter user, ID, condition
% 2. load in raw data
% 3. crop kTMP data to only include the 10 minutes during kTMP
% 4. create SDATA: 
%   a. antialiasing LPF filter (also exclude 60hz noise)
%   b. downsample
%   c. HPF
%   d. concatonate
% 5. inspect data and mark bad electrodes, if necessary
% 6. test various references
% 7. independent component analysis (ICA) to remove eye movements (credit to Dr. Christina Merrick and Dr. Assaf Breska)
% 8. semi-automatic artifact correction
% 9. interpolate bad electrodes, if necessary
% 10. save

%% 1. initialize 
user = 'owen'; %'owen' or 'christinamerrick'
ID = '07EL';
condition = 10; % kTMP AM freq
sham = true;

addpath(genpath(strcat('/Users/',user,'/Box/kTMP_EEG/BC-MR')));
addpath(genpath(strcat('/Users/',user,'/Box/EEG_pipeline')));
addpath(genpath(strcat('/Users/',user,'/Desktop/fieldtrip-master/fileio')));

if sham
    sham_str = '_sham';
else
    sham_str = '';
end

raw_data_path = strcat('/Users/',user,'/Box/kTMP_EEG/raw_data/',ID,'/',string(condition),'Hz_AM',sham_str); % The file selection dialog will open on this folder

% necessary functions:
% from fieldtrip: read_biosemi_bdf (plus all the read_24bit files), ReadBioSemiTriggerChannel
% from EEGlab (should have full toolbox as folder): runica, topoplot (plus '.locs' file, e.g. head71.locs)
% from LeonLab: HPF, LPF, multichanplot, remove_line_noise, rdi, get_neighbor_electrodes
%% 2. choose data file, load and append

% select files of a single subject
[data_file_list, ~] = uigetfile({'*.eeg','EEG file (*.eeg)' ;'*.*', 'All files (*.*)'},['Select raw file/s for S' num2str(ID)],raw_data_path,'multiselect','on');
% data_file_list = {'CM_test2.eeg'};
if isscalar(data_file_list)
    % No file was selected. 
else
    if ~iscell(data_file_list)        
        data_file_list = {data_file_list}; % Single file was selected - put it in a 1x1 cell array
    end
end

nSessions = length(data_file_list);

SDATA = struct;
SDATA.info = struct;
SDATA.metadata = struct;
SDATA.events = struct;

%% 3. Crop the kTMP .eeg file to the 10 minutes of stimulation (remove buffer time before/ after)

% load data
fileName = char(fullfile(raw_data_path,char(strcat(ID,"_kTMP.eeg"))));
kTMP_hdr = ft_read_header(fileName);
kTMP_session_data = ft_read_data(fileName);
kTMP_session_data = single(kTMP_session_data)';

%% mark the kTMP artifact when stimulation is turned on
kTMP_art_on = multichanplot(kTMP_session_data,'srate',kTMP_hdr.Fs,'channelnames', kTMP_hdr.label);

%% mark the 10 minutes of kTMP
if sham
    kTMP_mask = zeros(length(kTMP_session_data),1);
    ten_min_len = 10*60*kTMP_hdr.Fs;
    buffer = 10*kTMP_hdr.Fs;
    endsham = size(kTMP_session_data,1)-buffer;
    startsham = endsham-ten_min_len;
     
    kTMP_mask(startsham:endsham,:) = 1;
else
    kTMP_mask = multichanplot(kTMP_session_data,'srate',kTMP_hdr.Fs,'channelnames', kTMP_hdr.label);
end
%% mark the kTMP artifact when stimulation is turned off
kTMP_art_off = multichanplot(kTMP_session_data,'srate',kTMP_hdr.Fs,'channelnames', kTMP_hdr.label);

%% save kTMP artifacts
SDATA.kTMP_art = struct;

if exist('kTMP_art_on', 'var')
    SDATA.kTMP_art.on = kTMP_session_data(kTMP_art_on,:);
end

if exist('kTMP_art_off', 'var')
    SDATA.kTMP_art.off = kTMP_session_data(kTMP_art_off,:);
end

% crop and save kTMP data
kTMP_session_data = kTMP_session_data(kTMP_mask,:);
duration = size(kTMP_session_data, 1)/(kTMP_hdr.Fs * 60); % minutes
disp(strcat("kMTP duration: ", string(duration), " minutes")); 

%% 4. Create SDATA
EEGdata = [];
event_times = zeros(length(data_file_list), 2);
stop = 0;
data_times_check = [];
impedance_list = zeros(32,7);

% run over sessions
for sessIdx = 1:nSessions
    
    fileName = char(fullfile(raw_data_path,data_file_list{sessIdx}));

    if strcmp(data_file_list{sessIdx}, strcat(ID,"_kTMP_REAL.eeg"))
        hdr = kTMP_hdr;
        session_data = kTMP_session_data;
        disp("using cropped ktmp data (_kTMP_REAL)");    
    elseif strcmp(data_file_list{sessIdx}, strcat(ID,"_kTMP.eeg"))
        hdr = kTMP_hdr;
        session_data = kTMP_session_data;
        disp("using cropped ktmp data (_kTMP)");
    else
        hdr = ft_read_header(fileName);
        session_data = ft_read_data(fileName);
        session_data = single(session_data)';
    end
    
    disp(fileName)
    
    % a. Anti-aliasing low pass filter: The default filter is a 4th-degree non-causal Butterworth filter. 
    Fs = hdr.Fs;   % aquisition sampling rate
    lpf_freq = 50; % must be less than or equal to desired_sr/2
                   % lpf_freq < 60Hz will also exclude electrical noise.
    filt_deg = 4;  % default degree
    
    session_data = LPF(session_data,Fs,lpf_freq,filt_deg);
    
    % b. Down sample
    sr = 500; % desired sampling rate
    dsFactor = Fs/sr; % down sample factor
    session_data = downsample(session_data, dsFactor);
    
    % c. High pass filter: Exclude slow, non-neurological changes over changes over time (sweat, movement, ...).
    hpf_freq = 1;
    
    session_data = HPF(session_data,sr,hpf_freq,filt_deg);
    
    % d. concatonate
    EEGdata = [EEGdata; session_data];
    data_times_check = [data_times_check; size(session_data, 1)];
    impedance_list(:,sessIdx) = hdr.orig.impedances.channels;
    
    start = stop + 1;
    stop = stop + size(session_data, 1);
    event_times(sessIdx,:) = [start, stop];
    
    data_len = size(session_data, 1);
    disp(strcat("Duration: ", string(data_len/(sr*60)), " minutes"));
end

% log info
SDATA.data=EEGdata;
SDATA.events.event_times = event_times;
SDATA.events.event_times_labels = ["pre1","start","stop";"pre2","start","stop";"kTMP","start","stop";"post1","start","stop";"post2","start","stop";"post3","start","stop";"post4","start","stop"];

SDATA.info.id = ID;
SDATA.info.sampling_rate  = sr;
SDATA.info.channel_labels = hdr.label;
SDATA.info.impedancesAve  = mean(impedance_list,2);
SDATA.info.impedances     = impedance_list;

SDATA.info.channel_nums = SDATA.info.channel_labels;
for chan=1:length(SDATA.info.channel_nums)
    SDATA.info.channel_nums{chan}=[SDATA.info.channel_nums{chan} ' (' num2str(chan) ')'];
end

table_check = SDATA.events.event_times(:,2)-SDATA.events.event_times(:,1) + 1;
check = (table_check == data_times_check);  


%% 5. inspect data and mark bad electrodes

%% step 1: inspect data
multichanplot(SDATA.data, 50,'srate', sr, 'channelnames', SDATA.info.channel_nums);
disp((SDATA.events.event_times(3,2)-SDATA.events.event_times(3,1))/(SDATA.info.sampling_rate*60));

% check FFT of specific electrode
block = 1;
elec = 'F7';
x_limit = [0,70];
y_limit = [0,3];

start = SDATA.events.event_times(block,1);
stop  = SDATA.events.event_times(block,2);

elec_data = SDATA.data(start:stop,strcmp(SDATA.info.channel_labels,elec));
figure;
[spectrum,x_axis] = plotFFT(elec_data,SDATA.info.sampling_rate);
plot(x_axis,spectrum);
xlim(x_limit);
ylim(y_limit);
title(strcat(string(elec)," raw data"));

%% step 2: mark bad electrodes
bad_elec=input('input numbers of bad electrodes (in [], separated by commas): ');

if isfield(SDATA.metadata, 'bad_electrodes')
    SDATA.metadata.bad_electrodes=unique([SDATA.metadata.bad_electrodes, bad_elec]);
else
    SDATA.metadata.bad_electrodes=bad_elec;
end
SDATA.metadata.good_electrodes=setdiff(1:length(SDATA.info.channel_labels), SDATA.metadata.bad_electrodes);


%% 6. reference

%FCZ reference
FCZ_test = SDATA.data;
% multichanplot(FCZ_test,'srate', sr, 'channelnames', SDATA.info.channel_labels)

% common average reference
CAR_test = bsxfun(@minus,SDATA.data,mean(SDATA.data(:,SDATA.metadata.good_electrodes),2));
% multichanplot(CAR_test,'srate', sr, 'channelnames', SDATA.info.channel_labels)

% nose reference
nose_ref_test = bsxfun(@minus,SDATA.data(:,1:end),SDATA.data(:,strcmp(SDATA.info.channel_labels, 'ECG')));
% multichanplot(nose_ref_test,'srate', sr, 'channelnames', SDATA.info.channel_labels)

% reference elec of choice
elec = 'Oz';
elec_ref_test = bsxfun(@minus,SDATA.data,SDATA.data(:,strcmp(SDATA.info.channel_labels, elec)));
% multichanplot(elec_ref_test,'srate', sr, 'channelnames', SDATA.info.channel_labels)

%% save data before ICA, incase we want to go back and make changes to component selection
cd(strcat('/Users/owen/Box/kTMP_EEG/SDATA/',ID));
save(strcat(string(condition),'hz_SDATA_pre_ICA'),'SDATA')

SDATA_noseRef = SDATA;
SDATA_noseRef.data = nose_ref_test;
save(strcat(string(condition),'hz_SDATA_pre_ICA_nose_ref'),'SDATA_noseRef')
%% 7. ICA

%% step 1: prepare for ICA

manual_select_train_set=1;
train_set_range_secs=[10 130];

% remove bad electrodes
dataForICA=SDATA.data(:,SDATA.metadata.good_electrodes);
clean_channel_labels=SDATA.info.channel_labels(SDATA.metadata.good_electrodes);

%stage 3: define training set using multichanplot function
if ~manual_select_train_set
    ica_train_set=SDATA.data(train_set_range_secs(1)*sr:train_set_range_secs(2)*sr,:);
else
    
    show_data = [dataForICA];
    
    train_set_idx=multichanplot(show_data,50, 'srate', sr, 'channelnames', clean_channel_labels);
    ica_train_set=dataForICA(train_set_idx,:);
    
end

%% step 2: run ICA 

[weights,sphere,compvars,bias,signs,lrates,activations]  = runica(ica_train_set','pca',size(ica_train_set,2)-1);

unmix = weights*sphere;   % comps x chans
mix = pinv(unmix);

if ~isreal(mix)
    error('Warning: bad mixing matrix (complex values). Try reducing training data dimensions');
end

SDATA.ica.unmix = unmix;
SDATA.ica.mix = mix;

 %% step 3: inspect ICA component time course
C = (SDATA.ica.unmix * dataForICA')';

clc
close all

component_labels_with_bipolar=cell(size(C,2),1);
for comp=1:length(component_labels_with_bipolar)
    component_labels_with_bipolar{comp}=num2str(comp);
end

multichanplot(C, 50,'srate', sr);

clc

%% step 4: inspect ICA component topography
%enter ICA component to plot
compsToPlot = [1,15,29,11,13,14];
rows = 2;
numComps = length(compsToPlot);
figure;

for i = 1:numComps
    
    subplot(rows,3, i)
    compToPlot = compsToPlot(i);
    compTopo=zeros(length(SDATA.info.channel_labels),1);
    compTopo(SDATA.metadata.good_electrodes)=SDATA.ica.mix(:,compToPlot);
    compTopo(SDATA.metadata.bad_electrodes)=NaN;

    topoplot(compTopo,'headMR32.locs','electrodes','labels','style','map','shading','interp');
    colorbar;
    title(num2str(compToPlot));
end
%% step 5: remove components and reconstruct data
clc
remove_comp=input('input numbers of components to remove (in [], separated by commas): ');

Ncomps = size(SDATA.ica.mix,2);
cmp = true(Ncomps,1);
cmp(remove_comp) = false;
reconstructedData = ( SDATA.ica.mix * diag(cmp) * SDATA.ica.unmix * dataForICA' )';

%% step 6: Inspect reconstruction

% clc
% close all

multichanplot(SDATA.data, 50, 'srate', sr, 'channelnames', SDATA.info.channel_labels);
title('Pre-ICA Timeseries')

multichanplot(reconstructedData, 50, 'srate', sr, 'channelnames', SDATA.info.channel_labels(SDATA.metadata.good_electrodes));%, 'channelnames', SDATA.info.channel_labels);
title(strcat("Post-ICA Timeseries: ", string(remove_comp)));
%% step 7: approve reconstruction

acceptReconst=input('Accept reconstruction? (1=yes): ');
if acceptReconst==1
    SDATA.data(:,SDATA.metadata.good_electrodes) = single(reconstructedData);
end


%% 8. semi-automatic artifact rejection %%

%% step 1: automatic artifact detection
lpf_cutoff_for_rdi = 50;

abs_amp_threshold = 16384; % 16.384 mV = 16384 uV
low_act_allowed = 0.5; % amplitude threshold for low activity (electrode is flat for too long)
low_act_interval = 100; % ms
minmax_allowed = 500;
minmax_interval = 100; % 

sr = SDATA.info.sampling_rate;
art_margin = 0.2;
art_margin_sp=round(art_margin*sr/2);


channelsForRDI=SDATA.metadata.good_electrodes;                          
lpf_data_for_RDI = LPF(double(SDATA.data),sr,lpf_cutoff_for_rdi,3);

tic; amp_art = rdi(lpf_data_for_RDI,'channels',channelsForRDI,'m',[minmax_allowed minmax_interval art_margin_sp art_margin_sp], 'e',[-abs_amp_threshold abs_amp_threshold art_margin_sp art_margin_sp], 'l',[low_act_allowed low_act_interval art_margin_sp art_margin_sp]); toc
amp_art = any(amp_art,2);

%% step 2: inspect artifacts and modify marking

% easier to see artifact on non-preprocessed data, so load in raw data with the same timing to create mask
dispEEGdata = [];

% run over sessions
for sessIdx = 1:nSessions
    
    fileName = char(fullfile(raw_data_path,data_file_list{sessIdx}));

    if strcmp(data_file_list{sessIdx}, strcat(ID,"_kTMP_REAL.eeg"))
        hdr = kTMP_hdr;
        session_data = kTMP_session_data;
        disp("using cropped ktmp data (_kTMP_REAL)");    
    elseif strcmp(data_file_list{sessIdx}, strcat(ID,"_kTMP.eeg"))
        hdr = kTMP_hdr;
        session_data = kTMP_session_data;
        disp("using cropped ktmp data (_kTMP)");
    else
        hdr = ft_read_header(fileName);
        session_data = ft_read_data(fileName);
        session_data = single(session_data)';
    end
    
    disp(fileName)
    
    % save data
    dispEEGdata = [dispEEGdata; session_data];
end

channel_labels_RDI=SDATA.info.channel_labels(channelsForRDI);
amp_art_inspected = multichanplot(dispEEGdata(:,channelsForRDI),20,'srate',Fs,'markdata',amp_art,'channelnames',channel_labels_RDI);

% downsample the mask to match the SDATA
dsFactor = Fs/sr;
amp_art_inspected_test = downsample(amp_art_inspected, dsFactor);
multichanplot(SDATA.data(amp_art_inspected,channelsForRDI),20,'srate',sr,'markdata',amp_art,'channelnames',channel_labels_RDI);

%% step 3: accept artifact marking
acceptInspection=input('Accept artifact marking? (1=yes): ');
if acceptInspection==1
    SDATA.metadata.artifacts = amp_art_inspected;
end

%% save the data
% if not(exist(strcat('/Users/owen/Box/kTMP_EEG/SDATA/',ID),'dir'))
%     mkdir /Users/owen/Box/kTMP_EEG/SDATA
% end

cd(strcat('/Users/owen/Box/kTMP_EEG/SDATA/',ID));
SDATA.hdr = hdr;
SDATA.condition = condition;

save(strcat(ID,'_',string(condition),'hz','_SDATA'), 'SDATA', '-v7.3')

%% 9. Interpolate bad electrodes

% linear interpolation using lookup table
elecs = ["T7","T8"];

chans = SDATA.info.channel_labels;
neighborsIdx = zeros(size(chans,1),length(elecs));
elecsIdx = zeros(length(elecs),1);
neighbors = struct;

neighbors.Fp1 = ["F7","F3","Fz","Fp2"];
neighbors.Fp2 = ["F8","F4","Fz","Fp1"];
neighbors.F3 = ["Fz","Fc1","Fc5","F7","Fp1"];
neighbors.F4 = ["Fz","Fc2","Fc6","F8","Fp2"];
neighbors.T7 = ["F7","Fc5","C3","Cp5","Tp9"];
neighbors.T8 = ["F8","Fc6","C4","Cp6","Tp10"];
neighbors.Tp9 = ["T7","Cp5","P7"];
neighbors.Tp10 = ["T8","Cp6","P8"];

interpData = SDATA.data;

% get indicies for neighboring electrodes for interpolation
for i = 1:length(elecs)
    elec = elecs(i);
    disp(elec)
    elecsIdx(i,1) = find(strcmp(elec,SDATA.info.channel_labels));
    
    for j = 1:length(neighbors.(elec))
        neighbor = neighbors.(elec)(j);
        if ~ismember(neighbor,neighbors.(elec))
        neighborsIdx(:,i) =  neighborsIdx(:,i) + strcmp(chans,neighbor);  
        end
    end    
end

% update data to use interpolated electrodes
for i = 1:length(elecs)
    neighborsIdx = logical(neighborsIdx);
    interpData(:,elecsIdx(i,1)) = mean(SDATA.data(:,neighborsIdx(:,i)),2);
end

multichanplot(SDATA.data,20,'srate', sr, 'channelnames', SDATA.info.channel_labels)
multichanplot(interpData,20,'srate', sr, 'channelnames', SDATA.info.channel_labels)

%% approve interpolation
acceptInspection=input('Accept interpolation? (1=yes): ');

if acceptInspection==1
    SDATA.data = interpData;
    SDATA.interpolated_elecs = elecs;
end

%% 10. save the data
cd(strcat('/Users/owen/Box/kTMP_EEG/SDATA/',ID));
SDATA.hdr = hdr;
SDATA.condition = condition;

save(strcat(ID,'_',string(condition),'hz','_SDATA_interp'), 'SDATA', '-v7.3')