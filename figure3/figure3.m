load ../baseline/baseline.mat

%% set the colors: 
signalcolors = get(gca,'ColorOrder');
signals = {'control', 'CHIR', 'DLL1', 'PDGF', 'BMP4'};

test = brewermap(40,'Set2');
newcollist = [test(3,:); test(4,:); test(2,:)];

%% Merge data structures
% 
mergeraw = mergestructures(rawmult1,rawmult2);
mergelog= mergestructures(logmult1,logmult2);
mergefilt = mergestructures(filtmult1,filtmult2);
mergeH = mergestructures(Hmult1,Hmult2);

%% run the log normalization again with the correct samples
mergelog2 = lognorm(mergeraw);
mergefilt2 = filterdata(mergelog2,mergeraw,'method', 'hvg', 'thresh', hvgthresh, 'removeRibo', true);
mergeH2 = projectdata(feats, mergefilt2);

%% Build GMMs
rlist = 4;
sampsize = 1000;
mergegmms2 = buildgmm(mergeH2, sampsize, rlist, 1,1);

%% Plot heatmap of populations using fig1 genes
load ../Fig1-heatmap/proggenes.mat

gidx = progidx;
plotGMMheatmaplow(mergegmms2.gmmlist{37}.cluster(mergeH2.Hlist{37}'), mergelog2.logdflist{37}, mergelog.genes, true, gidx,false,[])
%print('figures/control_model_subpop_heatmap.pdf', '-dpdf', '-r300');

%% Fig SI Ranking. Generate figure to rank the data: 
% 
% keepinds = find(gmmdatareps.sampnums==6);
% keepinds = keepinds(1:5);
% vzd0gmm = gmmdatareps.gmmlist(keepinds);
% vzd0names =gmmdatareps.modelnames(keepinds);
% vzd0H = Hdata.Hlist(6);
% vzd0raw = rawdata.dflist(6);
% 
% newmergegmms = mergegmms;
% newmergegmms.gmmlist = [vzd0gmm, newmergegmms.gmmlist];
% newmergegmms.modelnames = [repmat({'SVZ-D0'},1,5),newmergegmms.modelnames];
% newmergegmms.dfnames = [{'SVZ-D0'},newmergegmms.dfnames];
% newmergegmms.sampnums = [repmat(1,1,5),newmergegmms.sampnums + 1];
% newmergegmms.bestlist = [2,newmergegmms.bestlist+5]
% 
% newmergeH = mergeH;
% newmergeH.Hlist = [vzd0H,mergeH.Hlist];
% newmergeH.dfnames = [{'SVZ-D0'},newmergeH.dfnames];
% newmergeH.sampnums = [1, newmergeH.sampnums+1];
% 
% newmeta = mergeraw.meta;
% newmeta{:,2} = newmeta{:,2} + 1; 
% newrow = [{'temp'}, num2cell([1,zeros(1,7)])];
% newmeta = [newmeta; newrow];

%%%%%%%%%%%%%% Add multiple samples

keepinds = [6,7,8,9,10,5]
vzd0gmm = gmmdatareps.gmmlist(keepinds);
vzd0names =gmmdatareps.modelnames(keepinds);
vzd0H = Hdata.Hlist(keepinds);
vzd0raw = rawdata.dflist(keepinds);

newmergegmms = mergegmms2;
newmergegmms.gmmlist = [vzd0gmm, newmergegmms.gmmlist];
newmergegmms.modelnames = [vzd0names,newmergegmms.modelnames];
newmergegmms.dfnames = [vzd0names,newmergegmms.dfnames];
newmergegmms.sampnums = [1:6,newmergegmms.sampnums + 6];
newmergegmms.bestlist = [1:6,newmergegmms.bestlist+6]

newmergeH = mergeH2;
newmergeH.Hlist = [vzd0H,mergeH.Hlist];
newmergeH.dfnames = [vzd0names,newmergeH.dfnames];
newmergeH.sampnums = [1:6, newmergeH.sampnums+6];

newmeta = mergeraw.meta;
newmeta{:,2} = newmeta{:,2} + 6; 
newrow1 = [{'temp'}, num2cell([1,zeros(1,7)])];
newrow2 = [{'temp1'}, num2cell([2,zeros(1,7)])];
newrow3 = [{'temp2'}, num2cell([3,zeros(1,7)])];
newrow4 = [{'temp3'}, num2cell([4,zeros(1,7)])];
newrow5 = [{'temp4'}, num2cell([5,zeros(1,7)])];
newrow6 = [{'temp5'}, num2cell([6,zeros(1,7)])];
newmeta = [newmeta; newrow1;newrow2;newrow3;newrow4;newrow5;newrow6];

[h,I] = rankGMMs(1, newmergegmms, newmergeH, newmeta, 100, 'heatmap', []);
set(gcf,'Position',[440   380   537   418])
print('figures/ranking_from_SVZ-d0.pdf', '-dpdf', '-r300');

%% Fig 3c: Plot model renderings
plotmodelPCA(newmergegmms, newmergeH, 1:2, 50, true)

%% Fig 3a: Plot the control renderings for all the data
% Generate a control model based on all the data

% picked specific samples in order to represent all populations evenly
allH = [mergeH.Hlist{[1,17,18,32]}]';
alllog = [mergelog.logdflist{[1,17,18,32]}];

N = min(5000,size(allH,1));
sampidx = randsample(size(allH,1),N);

sampH = full(allH(sampidx,:));
samplog = full(alllog(:,sampidx));
[coeff,score, LATENT, TSQUARED, EXPLAINED, pcamean]=pca(sampH);

% build a model from the sampled data: 
trainH = sampH(1:(round(N/2)),:)';
validH = sampH((round(N/2)+1):end,:)';
options = statset('MaxIter',300);

mkdir 'gmms/gmms_control/'
makegmms(trainH, validH, 'control', 'gmms/gmms_control/', 4, options, 1:4000, 4001:5000)
load 'gmms/gmms_control/control_4_gmm.mat'
controlgmm = gmfitall;

clusts = controlgmm.cluster(sampH);
figure; 
scatter3(-score(:,1),score(:,2),score(:,3),20,clusts,'filled')
colormap(brewermap(4,'Set2'))
colorbar()
set(gca,'view',[21.3000   24.4000])
print('figures/control_by_subpop.pdf', '-dpdf', '-r300');

%% Plot heatmap of populations
gidx = findgeneinds(logdata.genes, {'FABP7', 'LGALS1', 'SERPINE2','SOX2','ID3','NDRG2', 'APOE', 'ALDOC', 'ALDH1L1', 'GFAP','SPARC', 'SPARCL1', 'SYT4', 'NEUROD6', 'GAP43', 'LY6H','SOX11','MLLT11', 'TUBB3','SCG5', 'SOX4','CD24A','STMN3',  'CKS2','CENPA','CDC20','CCNB1','UBE2C','AURKA','KPNA2','CDK1'})
plotGMMheatmaplow(controlgmm.cluster(allH), alllog, mergelog.genes, true, gidx,false,[])
print('figures/control_model_subpop_heatmap.pdf', '-dpdf', '-r300');

%% Plot heatmap of populations using fig1 genes
load ../Fig1-heatmap/proggenes.mat

gidx = progidx;
plotGMMheatmaplow(controlgmm.cluster(allH), alllog, mergelog.genes, true, gidx,false,[])
print('figures/control_model_subpop_heatmap.pdf', '-dpdf', '-r300');

%% Specify a metadata ordering that will work for most cases

mr = [1,4,7,2,5,8,3,6,9,38,34,40,39,28,22,...
    1,13,10,2,14,11,3,15,12,33,...
    16,17,18,32,...
    3,...
    25,35,31,29,23,...
    26,27,30,29,20,...
    19,21,24,23,20]


%% Specify an ordering for training cases: 
tr = find(sum(uniquemeta{:,2:end}'>0)<3)
tr(tr == 37)=[]; % remove outlier case

%% Specify metadata as an X vector (normalized)
uniquemeta = unique(meta(:,2:end));
nx = uniquemeta{:,2:end}; 
X = nx ./ max(nx); % rescale metadata by maximums

%% Plot number of cells per condition
numcells = cellfun(@(x) size(x,2), newmergeH.Hlist(7:end))
figure; ax1=subplot(5,1,1:4);
plot(numcells(mr),'k','LineWidth',2);
set(gca,'XTick',[],'FontSize',16)
ylabel('numCells');
ax2 = subplot(5,1,5);
imagesc(X(mr,:)');
set(gca,'XTick',[],'FontSize',16)
linkaxes([ax1,ax2], 'x')
colormap(brewermap(20,'BuPu'))
set(gca, 'YTick', 1:length(ylabels), 'YTickLabel', ylabels, 'TickLabelInt', 'none')

%% SI Fig: Plot points in 3D for all models
nsamples = length(newmergeH.Hlist);
[r,c] = findsubplotsize(nsamples);

pclist = {};
figure; 
for i=1:nsamples
    currH = newmergeH.Hlist{i};
    currPC = (currH'-pcamean)*coeff;
    pclist{i} = currPC;
    currclust = controlgmm.cluster(currH');
    subplot(r,c,i);
    scatter3(currPC(:,1),currPC(:,2),currPC(:,3),10,currclust,'filled')
    colormap(brewermap(4,'Set2'))
    title(newmergeH.dfnames{i},'Interpreter', 'none')
    set(gca,'Xtick',[],'Ytick',[],'Ztick',[])
    set(gca,'FontSize',8)
end
set(gcf,'Position',[1466          84        1012         744])
print('figures/samples_scatter_by_subpop.pdf', '-dpdf', '-r300');

%% Query the data using sampled data: 

load 'gmms/gmms_control/control_4_gmm.mat'
controlgmm = gmfitall;

% designate the reference model  
refmodel = controlgmm;

numComponents = refmodel.NumComponents;
totaln = length(mergeH2.Hlist);
sampnums = mergeH2.sampnums;
cellnums = cellfun('size', mergeH2.Hlist, 2)

% construct a matrix of counts
allcounts = [];

for i=1:totaln % iterate over the samples

    currH = mergeH2.Hlist{i};
    % score the current dataset using the gmms:
    currPredict = refmodel.cluster(currH');

    sampcounts = [];
    for k=1:refmodel.NumComponents
        sampcounts(k) = sum(currPredict==k);
    end
    allcounts = [allcounts; sampcounts];
    
end

% rearrange counts in order to reflect pooling of pre-astrocytic
% progenitors
% 1 and 3 are AP
% 4 is CN --> 2
% 2 is NP --> 3

allcounts = [allcounts(:,1)+allcounts(:,3), allcounts(:,4),allcounts(:,2)]
means = (allcounts' ./ sum(allcounts'))'


%% BMP vs EGF for neuroblast and proliferative pre-astrocytic prog

popnames = {'AP','CN','NP'}
bmpidx=[9,6,3,8,5,2,7,4,1];

lowerlims = [0.1,0,0];
upperlims = [1,0.6,0.4];

for j=1:3
    bmpmeans = means(bmpidx,j);
    bmpmeans = reshape(bmpmeans,3,3);

    figure; 
    imagesc((bmpmeans))
    caxis([lowerlims(j),upperlims(j)])
    ylabel('BMP4')
    xlabel('EGF/FGF')
    title(popnames{j})
    colorbar('southoutside')
    set(gca,'XTick',1:3,'XTickLabel',[0.8,4,20])
    set(gca,'YTick',1:3,'YTickLabel',fliplr([0,4,20]))
    set(gcf,'pos',[ 440   570   162   228]);
    print(['figures/bmp_heatmap_',popnames{j},'.pdf'], '-dpdf', '-r300');
end

%% PDGF
pdgfidx = [1,10,13,2,11,14,3,12,15];
    
lowerlims = [0.1,0,0];
upperlims = [1,0.6,0.4];
for j=1:3
    pdgfmeans = means(pdgfidx,j);
    pdgfmeans = reshape(pdgfmeans,3,3);

    figure; 
    imagesc((pdgfmeans))
     caxis([lowerlims(j),upperlims(j)])
    ylabel('PDGF')
    xlabel('EGF/FGF')
    title(popnames{j})
    colorbar('southoutside')
    set(gca,'XTick',1:3,'XTickLabel',[0.8,4,20])
    set(gca,'YTick',1:3,'YTickLabel',fliplr([0,4,20]))
    set(gcf,'pos',[ 440   570   162   228]);
    print(['figures/pdgf_heatmap_',popnames{j},'.pdf'], '-dpdf', '-r300');
end


%% CHIR
chiridx = [1,16,2,17,3,18]
for j=1:3
    chirmeans = (means(chiridx,j));
    chirmeans = reshape(chirmeans,2,3); 

    figure; 
    imagesc(fliplr(chirmeans))
    ylabel('CHIR')
    xlabel('EGF')
    title(popnames{j})
    set(gca,'XTick',1:3,'XTickLabel',[0.8,4,20])
    set(gca,'YTick',1:2,'YTickLabel',fliplr([0,3]))
    set(gcf,'pos',[ 440   570   162   228]);
    caxis([lowerlims(j),upperlims(j)])
    colorbar('southoutside')
    print(['figures/chir_heatmap_',popnames{j},'.pdf'], '-dpdf', '-r300');
end


%% Other signals
restidx = [37,25,37,26,37,33,37,19]
% ctrl, fgf9, gm-csf, pdgf, ifng

for j=1:3
    restmeans = (means(restidx,j)');
    restmeans = flipud(reshape(restmeans,2,4));
    figure; 

    subplot(1,4,1)
    imagesc(restmeans(:,1))
    title(popnames{j})
    ylabel('FGF9')
    set(gca,'YTick',1:2,'YTickLabel',fliplr([0,20]))
    caxis([lowerlims(j),upperlims(j)])
    colorbar('southoutside')

    subplot(1,4,2)    
    imagesc(restmeans(:,2))
    title(popnames{j})
    ylabel('GM-CSF')
    set(gca,'YTick',1:2,'YTickLabel',fliplr([0,20]))
    caxis([lowerlims(j),upperlims(j)])
    colorbar('southoutside')

    subplot(1,4,3)
    imagesc(restmeans(:,3))
    title(popnames{j})
    ylabel('PDGF')
    set(gca,'YTick',1:2,'YTickLabel',fliplr([0,20]))
    caxis([lowerlims(j),upperlims(j)])
    colorbar('southoutside')
    
    subplot(1,4,4)    
    imagesc(restmeans(:,4))
    title(popnames{j})
    ylabel('IFNG')
    set(gca,'YTick',1:2,'YTickLabel',fliplr([0,20]))
    caxis([lowerlims(j),upperlims(j)])
    colorbar('southoutside')
    
    
    set(gcf,'pos',[ 440   570   162   228]);
    print(['figures/rest_heatmap_',popnames{j},'.pdf'], '-dpdf', '-r300');
end



% %% -----------------------------------------------------------
% % Save counts and X2 data as a csv file
% % ------------------------------------------------------------

csvwrite('counts_n.csv', counts) % total cell counts after gmm model query
writetable(uniquemeta, 'signals_conc.csv') % concentrations
csvwrite('signals_x.csv', X) % normalized concentrations
