function EEG_preproc(expType,subj,outDir,highPass,lowPass,causal,reref,epochStart,epochEnd,baseStart,baseEnd,varargin)
% function EEGFilterSri(highPass,lowPass,filtMeth,causal,reref,stimulus)
% This function is used to pre-process EEG data 
% Inputs:
%   expType - is the name of the experimental protocol you are analyzing e.g. vwfa
%   subj - a string containing the subject ID
%   outDir - output directory where pre-processed data should be saved
%   highPass - high pass filter frequency
%   lowPass - low pass filter frequency
%   causal - 1 if using a (non-linear minimum phase) casual filter or 0 for non-causal filter
%   reref - 1 if you want to compute average reference or 0 for no re-referencing
%   epochStart - in seconds how much before each event you'd like to start the epoch e.g. -.2
%   epoch end - in seconds how much after each event you'd like to end the epoch e.g. .7
%   baseStart - in ms the beginning of the baseline period e.g. -200
%   baseEnd - in ms the end of the baseline period e.g. -2
%   varargin - optional input
%       varargin{1} - stimulus you want to analyze e.g. 2nd stim for RA task vs 1st for MVPA
% 
% The function assumes all raw data files are saved as S<Subj #>_<expType>00<run #>.raw e.g.S776_vwfa001.raw    
% 
% e.g. EEG_preproc('vwfa','776','/home/Sri/data/preproc',.1,30,1,1,1)


data = [];
trialInfo = [];
badChan = [];
for i = 1:6
    dataPath = fileparts(which(['S' subj '_' expType '00' num2str(i) '.raw'])); 
    elecFile = which('Hydrocel_GSN_128_1.0_TRIM_mod.sfp');
    EEG = pop_readegi(fullfile(dataPath,['S' subj '_vwfa00' num2str(i) '.raw']));
    EEG = pop_chanedit(EEG, 'load',{elecFile 'filetype' 'autodetect'});

    % Performs filtering here with default filter order. Causal or non-causal can be specified
    EEG = pop_eegfiltnew(EEG, highPass,0,[],0,0,0,causal);
    EEG = pop_eegfiltnew(EEG, 0,lowPass,[],0,0,0,causal);

    % remove even events which correspond to second word
    EEG.event = EEG.event(word:2:length(EEG.event));
    for j = 1:length(EEG.event)
        EEG.event(j).urevent = j;
    end

    % use clean_rawdata (EEGlab pluggin may have to download) to identify bad electrodes
    originalEEG = EEG;
    EEG = clean_rawdata(EEG, [], -1, [], [], 20, []);
    badChannels = find(~EEG.etc.clean_channel_mask);
    badChan{i} = badChannels;

    %% alternatively can use built-in EEGlab function
        % [~, bC{1}] = pop_rejchan(EEG, 'elec',[1:128] ,'threshold',5,'norm','on','measure','prob','freqrange',[1 250]);
        % [~, bC{2}] = pop_rejchan(EEG, 'elec',[1:128] ,'threshold',5,'norm','on','measure','kurt','freqrange',[1 250]);
        % [~, bC{3}] = pop_rejchan(EEG, 'elec',[1:128] ,'threshold',5,'norm','on','measure','spec','freqrange',[1 250]);
        % badChan{i} = unique([bC{1}, bC{2}, bC{3}]);
    %%


    % reinterpolate bad electrodes
    EEG = pop_select(originalEEG, 'nochannel', badChannels);
    EEG = pop_interp(EEG,originalEEG.chanlocs,'spherical');
    
    %Epoch data and remove baseline
    EEG = pop_epoch( EEG, {'DIN'}, [epochStart epochEnd]);
    EEG = pop_rmbase( EEG, [baseStart    baseEnd]);

    %% Could Run ICA here
        % EEG = pop_runica(EEG, 'extended',1,'interupt','on');
    %%
    
    %Artifact detection using all channels and a threshold of -80 to 80
    [EEG, Indexes] = pop_eegthresh(EEG, 1, [1:128], -80, 80, epochStart, epochEnd, 0,0); % events are not removed here, but Indexes contains bad trials for the current run
    EEG2 = pop_eegthresh(EEG, 1, [1:128], -80, 80, epochStart, epochEnd, 0,1); 
    %% Potentiall use other criteria to remove trials
        % remove trials with linear trend
        % remove trials with abnormal spectra
    %%

    goodTrials = true(size(EEG.event));
    goodTrials(Indexes) = false; % variable is the length of the total number of trials. false values indicate an artifactual trial
    
    % Re-reference to average reference scheme using good trials only
    if reref
        EEG2.nbchan = EEG2.nbchan+1;
        EEG2.data(end+1,:,:) = zeros(1, EEG2.pnts,size(EEG2.data,3));
        EEG2.chanlocs(1,EEG2.nbchan).labels = 'initialReference';
        EEG2 = pop_reref(EEG2, []);
        EEG2 = pop_select( EEG2,'nochannel',{'initialReference'});
    end
    
    %replace non-average referenced data with average reference data
    EEG.data(:,:,goodTrials) = EEG2.data;

    %concatenate good trials and data across runs
    trialInfo = cat(2,trialInfo,goodTrials);
    data = cat(3,data,EEG.data);
end

if ~isempty(varargin)
    save(fullfile(outDir,sprintf(['Subject_%s_%.1f_%d_%d_%d_stimulus_%d.mat']),subj,highPass,lowPass,causal,reref,word),'data','trialInfo','badChan')
else
    save(fullfile(outDir,sprintf(['Subject_%s_%.1f_%d_%d_%d_.mat']),subj,highPass,lowPass,causal,reref,word),'data','trialInfo','badChan')
end
