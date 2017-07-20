function [accu,tBins2] = TDT_Sensor_MVPA(w,l,r,varargin)

subj = {'794','801','859','880','882','893','895','898','900'};
rad =r;
start = 145;
finish = 251;
window = w;
lag = l;
testNum = 10;
foldThresh = 50;
numPerm = 1000;

for s = 1:length(subj)
    foo = load(['Subject_' subj{s} '_' varargin{1} '.mat']);
%     foo2 = load([subj{s} '_Average_Reference.mat']);
    foo = foo.data;
    subjDat(s) = foo;
end
[tBins,tBins2] = createTBins(start,finish,window,lag);

elecs = spotlightCreater(foo,rad,'ROI');

accu = zeros(length(subj),size(tBins,1),length(elecs),numPerm);
clear foo
for s = 4:length(subj)
    foo = zeros(size(tBins,1),length(elecs),numPerm);
    for p = 1:numPerm
        cfg_sri = gatherData(subj{s},subjDat(s),1:3,4:6);
        if p>1
            cfg_sri = createPerm(cfg_sri,p);
        end
        if exist('rI')
            cfg_sri = createFolds(cfg_sri,foldThresh,testNum);
        else
            cfg_sri = createFolds(cfg_sri,foldThresh,testNum);
        end

        %%  Create Training and Testing Folds
        trainDesignMat = zeros(size(cfg_sri.dat,1),cfg_sri.numFolds);
        for f = 1:cfg_sri.numFolds
            trainDesignMat(cfg_sri.trainIdx{f},f) = 1;
        end
        testDesignMat = zeros(size(cfg_sri.dat,1),cfg_sri.numFolds);
        for f = 1:cfg_sri.numFolds
            testDesignMat(cfg_sri.testIdx{f},f) = 1;
        end
        %% Set Parameters
        cfg = setCfg(cfg_sri,cfg_sri.numFolds,trainDesignMat,testDesignMat);
        tic
        parfor t = 1:size(tBins,1)
            passedData.data = double((cfg_sri.dat(:,:,tBins(t,1):tBins(t,2))));
            passedData.dim = [NaN,NaN,NaN];
            [passedData,cfg2] = fill_passed_data(passedData,cfg,cfg_sri.allLabels_overall.labels_supracat_clean);
            cfg2.results.output = {'accuracy'};
            foo(t,:,p) = decoding_sri(cfg2,passedData,elecs);
            close all
        end
        clc
        toc
    end
    save(['/home/sri/Documents/MATLAB/Scripts/Results/Sensor/perm_subj_' subj{s} '.mat'],'foo')
    accu(s,:,:,:) = foo;
    clear foo
end
