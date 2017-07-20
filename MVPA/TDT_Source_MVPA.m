function [accu,tBins2] = TDT_Source_MVPA(w,l,r,varargin)

subj = {'776','794','801','859','880','882','893','895','898','900'};
rad =r;
start = 156;
finish = 251;
window = w;
lag = l;
testNum = 10;
foldThresh = 100;
numPerm = 1000;
load('Kernel_232.mat')
subjKernel = fixKernel();
grid = InvKernel(1).GridLoc;
gridLoc.elec.pnt = grid;

for s = 1:1
    foo = load(['Subject_' subj{s} '_' varargin{1} '.mat']);
    foo = foo.data;
    subjDat(s) = foo;
end
[tBins,tBins2] = createTBins(start,finish,window,lag);
elecs = spotlightCreater(gridLoc,rad,'ROI');
tBins = [187 231];
tBins2 = [tBins2(15,2) tBins2(36,2)];
accu = zeros(length(subj),size(tBins,1),length(elecs),numPerm);
clear foo


for s = 1:1
    k = subjKernel{s};
    foo = zeros(size(tBins,1),length(elecs),numPerm);
    for p = 1:numPerm
        cfg_sri = gatherData(subj{s},subjDat(s),1:3,4:6,k);
        if p>1
            cfg_sri = createPerm(cfg_sri,p);
        end
        cfg_sri = createFolds(cfg_sri,foldThresh,testNum);

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
        t=1;
%         for t = 1:size(tBins,1)
            passedData.data = double(mean(cfg_sri.dat(:,:,tBins(t,1):tBins(t,2)),3));
            passedData.dim = [NaN,NaN,NaN];
            [passedData,cfg2] = fill_passed_data(passedData,cfg,cfg_sri.allLabels_overall.labels_supracat_clean);
            cfg2.results.output = {'accuracy'};
            foo(1,:,p) = decoding_sri(cfg2,passedData,elecs);
            close all
%         end
        clc
        toc
    end
    save(['/home/malone/Documents/MATLAB/Sri/Results/Source/perm_subj_' subj{s} '.mat'],'foo')
    accu(s,:,:,:) = foo;
end
