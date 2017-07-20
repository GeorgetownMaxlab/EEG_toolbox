function [coh_cond1,coh_cond2,tBins,tBins2] = runSeedCoh_Cond(subj,dataDir,filestem,outDir,freqoi,seed1,seed2,t1,t2,lag,numCycles,numPerm)
%function [coh_cond1,coh_cond2,tBins,tBins2] = runSeedCoh_Cond(subj,dataDir,filestem,outDir,freqoi,seed1,seed2,t1_idx,t2_idx,lag,numCycles,numPerm)
%This function computes coherence between two seeds at a particular frequency or frequency band in a series of sliding time windows
%It can also provide permuted coherence for later permutation and cluster-correction techniques
% Input:
%   - subj: cell array of subject IDs to be analyzed
%   - dataDir: location of the data
%   - filestem: suffix following subject ID e.g. for 776_0.2_30_sep_1_0_word_1.mat the filestem would be '_0.2_30_sep_1_0_word_1.mat'
%   - freqoi: cell array containing frequencies of interest. Each cell may contain a single frequency or an array of frequencies (for band coherence)
%   - seed1: array of channel indices to be used as seed e.g. 1:128 to calculate coherence to every sensor in the brain
%   - seed2: array of channel indices to be used as seed
%   - t1: start time in ms. Must be event for 500hz sampled data.
%   - t2: end time in ms. Must be event for 500hz sampled data.
%   - lag: lag in ms of sliding window smallest value can be 2ms for 500hz sampled data
%   - numCycles: number of cycles to be used to calculate coherence. Can use 1, but 2-3 is recommended
%   - numPerm: number of permutations you'd like to find. If no permutations are needed set equal to 1
%
%

% Read in all subject data
for s = 1:length(subj)
    foo = load(fullfile(dataDir,[subj{s} filestem]));
    foo = foo.data;
    Data(s) = foo;
    clear foo
end

% need to convert times into corresponding data indices
t1_idx = t1/2+101; 
t2_idx = t2/2+101;
lag = lag/2;

% iterate through frequencies
for f = 1:length(freqoi) 
    %calculate time required for a single cycle given a certain frequency
    freq = freqoi{f};
    if numel(freq)>1 % in case user wants to calculate coherence over a band
        timwinsample = getWinLength(min(freq),numCycles); % use the smallest frequency in the band
    else
        timwinsample = getWinLength(freq,numCycles);
    end
    % get time windows
    [tBins{f}, tBins2{f}] = getTimeWindow(t1_idx,t2_idx,timwinsample,lag);
    times = tBins{f};

    coh_1 = single(zeros(length(subj),length(seed1),length(seed2),size(times,1),numPerm));
    coh_2 = single(zeros(length(subj),length(seed1),length(seed2),size(times,1),numPerm));
    for s = 1:length(subj)
        tic
        % Format subject data into appropriate form
        [word_dat,allLabels] = data_preproc(Data(s),subj{s});
        % These lines only work for sri to filter data by condition. If we put data in single 
        % format perhaps an all-purpose function could be written to filter different conditions
        labels_1 = allLabels.labels_supracat_clean==1;
        labels_2 =allLabels.labels_supracat_clean==0;
        dat_1 = word_dat(labels_1,:,t1_idx:t2_idx);
        dat_2 = word_dat(labels_2,:,t1_idx:t2_idx);
        times = times-t1_idx+1; % have to make sure indices are consistent after selecting a subset of time points
        coh_1_tmp = zeros(length(seed1),length(seed2),size(times,1),numPerm);
        coh_2_tmp = zeros(length(seed1),length(seed2),size(times,1),numPerm);
        % Calculate coherence in each time window
        parfor t = 1:size(times,1)
            win = times(t,1):times(t,2);
            % loop through permutations 
            for p = 1:numPerm
                if p >1 % first index of output is actual result and is calculated with un-permuted data
                    [dat_tmp,dat_tmp]= permuteData(dat_1(:,:,win),dat_2(:,:,win));
                else
                    dat_tmp = dat_1(:,:,win);
                    dat_tmp = dat_2(:,:,win);
                end
                % Compute Cohrence
                coh_2_tmp(:,:,t,p) = computeCoherence(dat_tmp,seed1,seed2,win,freq,timwinsample);
                coh_1_tmp(:,:,t,p) = computeCoherence(dat_tmp,seed1,seed2,win,freq,timwinsample);
            end
        end
        toc
        coh_2(s,:,:,:,:) = coh_2_tmp;
        coh_1(s,:,:,:,:) = coh_1_tmp;
    end
    coh_cond1{f} = coh_1;
    coh_cond2{f} = coh_2;
end
end


function [dat_all,allLabels] = data_preproc(data,subject)
% Put data into matrix form
allLabels = labelAssign(data,subject,'both');
data2.trial = data.trial(allLabels.art_trails);
dat_all = single(zeros(length(data2.trial),size(data2.trial{1},1),size(data2.trial{1},2)));
for t =1:length(data2.trial)
    dat_all(t,:,:) = data2.trial{t}; % dat is your data matrix
end
clear t data data2 conds_Idx
end

function timwinsample = getWinLength(freq,numCycles)
fsample =500; % sampling rate
timwin = numCycles./freq;
timwinsample = round(timwin .* fsample);
end

function [tBins, tBins2] = getTimeWindow(start,finish,window,lag)
clear tBins clear tBins2
count = 1;
for t = start:lag:finish-window
    tBins(count,:) = [t t+window];
    count = count+1;
end
tBins2 = (tBins-101)*2;
clear count t
end

function [newData1,newData2]= permuteData(data1,data2)
combData = cat(1,data1,data2);
randIdx = randperm(size(combData,1));
newData1 = combData(randIdx(1:size(data1,1)),:,:);
newData2 = combData(randIdx(size(data1,1)+1:end),:,:);
end

function coh = computeCoherence(data,chan1,chan2,win,freq,timwinsample)
numTrials = size(data,1);
X = squeeze(data(:,chan1,:));
Y = squeeze(data(:,chan2,:));
X = permute(X,[3,1,2]);
Y = permute(Y,[3,1,2]);
X2 = reshape(X,size(X,1)*size(X,2),length(chan1));
Y2 = reshape(Y,size(Y,1)*size(Y,2),length(chan2));
X2 = repmat(X2,1,length(chan2));
Y2 = repmat(Y2,1,length(chan1));
idx = repmat([1:length(chan2)],1,length(chan1));
idx2 = repmat([1:length(chan1)],1,length(chan2));
[~,idx3] = sort(idx,'ascend');
Y2 = Y2(:,idx3);
Cxy= mscohere(X2,Y2,numTrials,[],freq,500);
Cxy = mean(Cxy);
coh = reshape(Cxy,length(chan1),length(chan2)).';
coh = coh';
end

