% This MATLAB script can be used to reproduce the result of n-back ERP
%
% Please download BBCItoolbox to 'MyToolboxDir'
% Please download dataset to 'NirsMyDataDir' and 'EegMyDataDir'
% The authors would be grateful if published reports of research using this code
% (or a modified version, maintaining a significant portion of the original code) would cite the following article:
% Shin et al. "Simultaneous acquisition of EEG and NIRS during cognitive tasks for an open access dataset",
% Scientific data (2017), under review.
% NOTE: Figure may be different from that shown in Shin et al. (2017) because EOG-rejection is not performed.

clear all; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%% modify directory paths properly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MyToolboxDir = fullfile('C:','Users','shin','Documents','MATLAB','bbci_public-master');
% % % MyToolboxDir = fullfile('C:','','shin','Documents','MATLAB','bbci_toolbox_latest_ver');

% % %WorkingDir = fullfile('C:','Users','shin','Documents','MATLAB','scientific_data');
% % % NirsMyDataDir = fullfile('F:','scientific_data_publish','rawdata','NIRS');

% EegMyDataDir = fullfile('F:','scientific_data_publish','rawdata','EEG','without EOG');
% % %EegMyDataDir = fullfile('F:','scientific_data_publish','rawdata','EEG');
MyToolboxDir = fullfile('G:','workspace','bbci_public');
WorkingDir = fullfile('G:','workspace','scientific_data');
EegMyDataDir = fullfile('G:','workspace','Capstone_490','Data');

StatisticDir = fullfile('G:','workspace', 'scientific_data', 'behavior');
NbackDir = fullfile(StatisticDir,'n-back','summary');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(MyToolboxDir);
startup_bbci_toolbox('DataDir',EegMyDataDir,'TmpDir','/tmp/','History',0); %did EegMyDataDir instead of NirsMyDataDir
cd(WorkingDir);

addpath(genpath(pwd));
addpath(UtilDir);

%% initial parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subdir_list = {'VP001-EEG'};
nback = 3; % 0 / 2 / 3 -back
ival_epo = [-0.1 1] * 1000;
ival_base = [-0.1 0] * 1000;
ival_scalps = [0 0.2; 0.2 0.4; 0.4 0.6; 0.6 0.8; 0.8 1]*1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for vp = 1 : length(subdir_list)
     disp([subdir_list{vp}, ' was started']);

    % correct epoch selection
    vpDir = fullfile(NbackDir, subdir_list{vp});
    cd(vpDir);
    load summary1; load summary2; load summary3;   
    
    summary = [summary1.result; summary2.result; summary3.result];
    summary = reshape(summary', 1, size(summary,1)*size(summary,2))';
    
    correctIdx = find(summary == 1);
    incorrectIdx = find(summary ~= 1);
     
    cd(WorkingDir);
 
    %% Load EOG-free data
    loadDir = fullfile(EegMyDataDir,subdir_list{vp});
    cd(loadDir);
    load cnt_nback; load mrk_nback; load mnt_nback;
    cd(WorkingDir);

    disp(['eeg was loaded dawg']);
    %% Select EEG channels only
    cnt_nback = proc_selectChannels(cnt_nback, 'not', '*EOG'); % remove EOG channels (VEOG, HEOG)    
    mnt_nback = mnt_setElectrodePositions(cnt_nback.clab);
    
    disp(['eog shawties get denied!']);
    %% 50 Hz mains filtering
    [b, a] = butter(3, [1 45]/cnt_nback.fs*2, 'bandpass'); % 1-45 Hz band-pass filtering
    cnt_nback = proc_filtfilt(cnt_nback, b, a);
    
    disp(['we filterin dem dawgs!']);
    %% Artifact rejection based on variance criterion
    mrk_nback = reject_varEventsAndChannels(cnt_nback, mrk_nback, ival_epo, 'verbose', 1);
    
    
   disp('artifact rejection tings!']);
    %% Segmentation
    epo = proc_segmentation(cnt_nback, mrk_nback, ival_epo);
    
    %% Select epoch with correct answer
    epo = proc_selectEpochs(epo, 'not', incorrectIdx);
    epo.xUnit = 'ms'; 
    epo.yUnit = '\muV';
    disp([subdir_list{vp},': ',num2str(length(incorrectIdx)), ' epoch(s) was/were rejected due to incorrect answer']);
    
    disp(['select appropriate epoch donezo!']);
    %% Select class for 0-back / 2-back / 3-back
    switch nback
        case 0
            epo = proc_selectClasses(epo, '0-back target');
        case 2
            epo = proc_selectClasses(epo, {'2-back non-target','2-back target'});
        case 3
            epo = proc_selectClasses(epo, {'3-back non-target','3-back target'});
    end
    
    if vp == 1
        epo_all = epo;
    else
        epo_all = proc_appendEpochs(epo_all, epo);
    end
    
    clear epo
    disp(['categorized em styllll!']);
end

%% Baseline correction
epo_all = proc_baseline(epo_all, ival_base);
disp(['fam the baseline is a go-go!']);
%% Calculate the pointwise difference between class means
epo_diff = proc_classmeanDiff(epo_all, 'Stats',1); 
epo_diff = proc_baseline(epo_diff, [-100 0]);

% Merge epo_all & epo_diff
epo_all = proc_appendEpochs(epo_all, epo_diff);

% Dimentionality correction
epo_all.t = epo_all.t(:,:,1);

disp(['calculated pointwise difference buddy!']);
%% Trial-Average
epo_all = proc_average(epo_all, 'Stats',1);

disp(['trial avg is set bro!']);
%% Display
% epo_all = proc_movingAverage(epo_all, 100); % for smoothing
epo_all.className = {'non-target','target','non-target - target'}; % add correct class lables

clab = {'Cz','Pz'};
%%%% THIS IS THE PART I CANNOT FIX %%%%
fig_set(3, 'Toolsoff', 0, 'Resize', [1 2]);
 H = plot_scalpEvolutionPlusChannel(epo_all, mnt_nback, clab, ival_scalps, defopt_scalp_erp, ...
    'Extrapolation', 0, 'Contour', 0, 'printival', 1, 'Colormap', cmap_posneg(101), 'PlotStat', 'sem', ...
    'yLim', [-1 1]*6, 'CLim', [-1 1]*8);