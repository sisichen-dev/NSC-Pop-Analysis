function gmmdata = buildgmm(Hdata, sampsize, rlist, iterations, refindex)
% gmmdata = buildgmm(Hdata, sampsize, rlist, iterations)
% returns structure with log-normalized data
%  gmmdata.expname = experiment name
%  gmmdata.dfnames = names of all samples
%  gmmdata.refname = name of the reference sample (pick first in list)
%  gmmdata.gmmlist = list of gmm models
k = size(Hdata.Hlist{1},1);
gmmfolder = ['gmms/gmms_',num2str(k),'/'];
Hlist = Hdata.Hlist;
dfnames = Hdata.dfnames;
refname = Hdata.refname;

fprintf('\n')
disp('************************************')
disp('Building gmm models for data....')
disp('************************************')
fprintf('\n')

% Check if gmmfolder already exists, if not make it.
if (exist(gmmfolder,'dir')==7)
    
else
    mkdir(pwd,gmmfolder)
end


%% Build gmms for all data
    
gmmlist={};
gmmnames={};
xBIClist = [];
unusedCellsList = {};
ctrainlist={};
cvalidlist={};

% only used to load tabula muris data
% xBIC=1
% trainIdx = [1:10]
% validIdx = [11:20]

for iter=1:iterations
    for i=1:length(dfnames)
    
        name = dfnames{i};
        sampfolder = [gmmfolder, num2str(i),'_',name,'/'];
        curritergmm = [sampfolder, name, '_iter', num2str(iter), '_gmm.mat'];
% only used to load tabula muris data
%         gmmfolder = 'gmms/'
%         sampfolder = [gmmfolder, name, '/']
%         curritergmm = [sampfolder, name, '_iter', num2str(iter), '_gmm.mat']
        
        if (exist(curritergmm,'file')==2)
            disp(['Iteration ', num2str(iter),' GMM already learned for : ', name]);
            load(curritergmm)
            gmmlist{end+1} = gmfitall;
            gmmnames{end+1} = name;
            xBIClist(end+1) = xBIC;
            ctrainlist{end+1} = trainIdx;
            cvalidlist{end+1} = validIdx;
        else
            mkdir(pwd,sampfolder)
            disp(['Making iteration ', num2str(iter),' gmm for sample : ', name]);
            
            testH = Hlist{i};

            % Downsample testH to prevent model overfitting:
            
            % If there are fewer cells than 2x sample size, then use a
            % smaller sample size
            if size(testH,2) < (2 * sampsize)
                s = size(testH,2) * 0.5;
                trainIdx = randsample(size(testH,2), s);
                validIdx = setdiff(1:size(testH,2), trainIdx);
            else % get 
                s = sampsize;
                trainIdx = randsample(size(testH,2), s);
                restIdx = setdiff(1:size(testH,2), trainIdx);
                validIdx = restIdx(randsample(length(restIdx), sampsize));
            end
            
            trainH = testH(:,trainIdx);
            validH = testH(:,validIdx);
            
            % make gmms for all 
            makegmms(trainH, validH, name, sampfolder, rlist, trainIdx, validIdx)
            
            % Load in best gmm based on BIC:
            [bestgmm, xBIC] = loadbestgmm(name,sampfolder, gmmfolder, rlist, iter, testH, trainIdx, validIdx);
            gmmlist{end+1} = bestgmm;
            gmmnames{end+1} = name;
            xBIClist(end+1) = xBIC;
            ctrainlist{end+1} = trainIdx;
            cvalidlist{end+1} = validIdx;

        end
    end
end


%% Compile gmmdata structure

gmmdata.expname = Hdata.expname;
%gmmdata.refname = anchorrefname;
%gmmdata.dfnames = reshape(repmat(dfnames, iterations,1),1,iterations*length(dfnames));
gmmdata.modelnames = gmmnames;
gmmdata.dfnames = dfnames;
gmmdata.gmmlist = gmmlist;
gmmdata.ctrainlist = ctrainlist;
gmmdata.cvalidlist = cvalidlist;
gmmdata.xBIClist = xBIClist;
gmmdata.sampnums = repmat(Hdata.sampnums, 1, iterations);

% Designate a single reference model from the controls
gmmdata.refindex = refindex;
gmmdata.bestlist = 1:length(dfnames);

gmmdata.refname = Hdata.refname;

% cluster the data and store within the structure
gmmdata = gmmcluster(gmmdata,Hdata);
% calculate subpopulation entropies and store within the structure
gmmdata = gmmEntropy(gmmdata);

function [bestgmm, xBIC] = loadbestgmm(name, sampfolder, gmmfolder, rlist, iter, testH,trainIdx, validIdx)

xBIClist = [];
BIClist = [];
for i = 1:length(rlist)
    r = rlist(i);
    load([sampfolder, name,'_', num2str(r), '_gmm.mat'])
    % BIClist(i)= gmfitall.BIC;
    % check condition number of matrices in gmfitall
    if checkcond(gmfitall.Sigma)
        BIClist(i) = calcBIC(gmfitall, testH(:,trainIdx));
        xBIC = calcBIC(gmfitall, testH(:,validIdx));
        xBIClist(i)= xBIC;
    else
        xBIC = Inf;
        BIClist(i) = Inf;
        xBIClist(i) = Inf;
    end    
end

% find index that has best model
ind = find(xBIClist==min(xBIClist));
disp(['The best gmm model for ', name, ' has ', num2str(rlist(ind)), ' components'])
load([sampfolder, name,'_', num2str(rlist(ind)), '_gmm.mat']);
bestgmm = gmfitall;

figure;
plot(rlist,BIClist,'b'); hold on
plot(rlist,xBIClist,'r'); hold on;
scatter(rlist(ind), xBIClist(ind), '*r'); hold on;
title(name);
print([gmmfolder, num2str(i),'_',name, '_modelselection.pdf'],'-dpdf','-r300')

% find index that has best model
disp(['The best gmm model for ', name, ' has ', num2str(rlist(ind)), ' components'])
load([sampfolder, name,'_', num2str(rlist(ind)), '_gmm.mat']);
save([sampfolder, name, '_iter', num2str(iter), '_gmm.mat'], 'gmfitall', 'xBIC','trainIdx', 'validIdx')
bestgmm = gmfitall;

close all;

function bool = checkcond(Sigma)
% function bool = checkcond(gmfitall.Sigma)
% This function checks the condition number to make sure it is relatively
% reasonable in size. Returns true if condition numbers are okay (small
% enough)
conds = [];
for i=1:size(Sigma,3);
    currcond = cond(Sigma(:,:,i));
    conds = [conds;currcond];
end
bool = ~sum(conds>1e10);