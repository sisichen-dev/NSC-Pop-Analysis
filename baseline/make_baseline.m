dfnames = {'E15','E18-9k','P4', 'AdCo', 'NSC2', 'SVZ-D0', 'SVZ-D3', 'SVZ-D5', 'SVZ-D7', 'SVZ-D18', 'SVZ-D5-CHIR', 'SVZ-D5-DLL', 'SVZ-D5-PDGF','SVZ-D5-BMP4'}
rawdata = loadraw('../data/', dfnames, refname,expname,'genes.tsv','mtx');

load rawdata.mat
logdata = lognorm(rawdata);
genes = logdata.genes;

%% Normalize all data
logdata = lognorm(rawdata);

%% Remove RBCs: 
[logdata, rawdata] = removeRBCs(logdata, rawdata, 'mouse');

%% Filter the data
hvgthresh = 1.4;
filtdata = filterdata(logdata, rawdata, 'method', 'hvg', 'thresh', hvgthresh, 'removeRibo', true);

% %% Check what the PC data looks like: 
% pcafeats = makePCAfeat(filtdata,5000,10);
% plotPC3D(pcafeats, collist);

%% Now run standard popAlign algorithm: 
klist = 12;
sampsize = 100;
feats = makefeat(filtdata, klist, sampsize, 'features/');

%% Project data 
Hdata = projectdata(feats, filtdata);

%% Build gmms
rlist = 1:10;
gmmdata = buildgmm(Hdata, 1000, rlist, 1);
gmmdatareps = buildgmm(Hdata, 5000, rlist, 20);

plotmodelPCA(gmmdatareps,Hdata,1:2,100,false);

%% Calculate tsne
Hdata = calcTSNE(Hdata);

%% Check that the gmms here look good: 
Hdata = checkGMMs(gmmdatareps,Hdata, logdata,feats);

%%
%-------------------------------------
% Load in the screens
%-------------------------------------

rawmult1 = loadscreen('../data/', 'NPC-MULT-1', 'EGF_FGF-20', 'mult-1','genes.tsv','mtx');
rawmult2 = loadscreen('../data/', 'NPC-MULT-2', 'EGF_FGF-20', 'mult-2','genes.tsv','mtx');

%% Remove all data with less than minN number of cells: 
minN = 90;
rawmult1 = cullsamples(rawmult1,minN);
rawmult2 = cullsamples(rawmult2,minN);

%% Normalize all data
logmult1= lognorm(rawmult1);
logmult2= lognorm(rawmult2);

%% Find features
filtmult1 = filterdata(logmult1, rawmult1, 'method', 'hvg', 'thresh', hvgthresh) ; 
filtmult2 = filterdata(logmult2, rawmult2, 'method', 'hvg', 'thresh', hvgthresh) ; 

%% Project data
Hmult1 = projectdata(feats, filtmult1);
Hmult2 = projectdata(feats, filtmult2);

%% Remove methanol fixed sample from the multiplexed dataset:

rawmult1.dflist = rawmult1.dflist(1:18);
rawmult1.dfnames = rawmult1.dfnames(1:18);
rawmult1.sampnums = rawmult1.sampnums(1:18);

logmult1.dfnames =logmult1.dfnames(1:18);
logmult1.logdflist = logmult1.logdflist(1:18);
logmult1.sampnums = logmult1.sampnums(1:18);

filtmult1.filtdflist = filtmult1.filtdflist(1:18);
filtmult1.sampnums = filtmult1.sampnums(1:18);
filtmult1.dfnames = filtmult1.dfnames(1:18);

Hmult1.dfnames = Hmult1.dfnames(1:18);
Hmult1.Hlist = Hmult1.Hlist(1:18);
Hmult1.sampnums = Hmult1.sampnums(1:18);

%% Save all data in a baseline file. 
save('baseline.mat','-v7.3')

%-
