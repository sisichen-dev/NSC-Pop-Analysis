
%% load resource data
load ../data/proggenes.mat

%% load only the E18, P4, and SVZ-D0 datasets and make PCA features based on this data

raw1 = loadraw('../data/', {'AdCo','P4','E18','SVZ-D0'}, 'SVZ-D0','init','genes.tsv','mtx');
log1 = lognorm(raw1);
[log1, raw1] = removeRBCs(log1, raw1, 'mouse'); % remove RBCs
hvgthresh = 1.6;
filt1 = filterdata(log1, raw1, 'method', 'hvg', 'thresh', hvgthresh, 'removeRibo', true);
klist = 10;
sampsize = 500;
[pcafeats, pcaH1] = makePCAfeat(filt1, sampsize, 10);

%% Batch correct the data using coherent point drift
dflist = cellfun(@(x) x', pcaH1.Hlist, 'Uniform', false);
[newdflist, transforms] = daisy_CPD(dflist,500, 1,{'AdCo','P4','E18','SVZ-D0'});

%% Make new pcaH1
pcaH1_bc = pcaH1;
pcaH1_bc.Hlist = cellfun(@(x) x',newdflist,'Unif',false);

%% Build gmms
rlist = 1:9;
gmmdata = buildgmm(pcaH1_bc, 1000, 1:5, 1,1);

%% -------------------------------------------------------------
% Fig 1a: Scatter plot of stages
% -------------------------------------------------------------

figure;
collist = viridis(3);
% combine all data and make a color vector
for i=1:3
    currpcs = pcaH1_bc.Hlist{i};
    scatter(currpcs(1,:),currpcs(2,:),8,collist(i,:),'filled','MarkerFaceAlpha',0.5);hold on;

end
plotinds = randsample(1:size(all_c,2),size(all_c,2),false);

collist = viridis(3);
% combine all data and make a color vector
allpcs = [pcaH1_bc.Hlist{:}];
all_c =[]
for i=1:3
    all_c = [all_c, repmat(i,1,length(pcaH1_bc.Hlist{i}))];hold on;
end
plotinds = randsample(1:size(all_c,2),size(all_c,2),false);

% scatter points randomly so that the different stages overlay 
figure; 
scatter(allpcs(1,plotinds),allpcs(2,plotinds),8,all_c(plotinds),'filled','MarkerFaceAlpha',0.5);hold on;
colormap(collist);
% xlabel('PC1')
% ylabel('PC2')
set(gca,'XTick',[0,5,10,15])
set(gca,'YTick',[0,5,10,15])
set(gca,'FontSize',16)
axis tight;
set(gcf,'pos',[360   519   310   262]);
print('figures/Fig1_stages.png','-dpng','-r300');

%% -------------------------------------------------------------
% Fig 1a inset: Scatter plot by genes
% -------------------------------------------------------------
diffgenes = {'SYT1','SOX10','ALDOC'}; % OR GFAP
gene = 'ALDOC';
gidx = find(strcmp(log1.genes,gene));

figure;
for j=1:length(diffgenes)
    subplot(1,length(diffgenes),j);
    gene = diffgenes(j);
    gidx = find(strcmp(log1.genes,gene));
    for i=1:3
        currpcs = newdflist{i}';
        currcol = log1.logdflist{i}(gidx,:);
        scatter(currpcs(1,:),currpcs(2,:),1,sqrt(currcol),'filled','MarkerFaceAlpha',1); hold on;
        colormap(magma);
    end
    colorbar
    title(gene)
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    axis tight;
end
set(gcf,'pos',[122   602   636   125]);
print('figures/Fig1_diffgenes.pdf', '-dpdf', '-r300')

%% -------------------------------------------------------------
% Fig 1b: Plot model PCA 
% -------------------------------------------------------------
plotmodelPCA(gmmdata, pcaH1, 1:2, 100,'PCA', true, 'weighted', true)


%% -------------------------------------------------------------
% Fig 1c: Scatter plot by subpopulation
% -------------------------------------------------------------
collist = [128,128,128;
    226,75,51;
    59, 82,  139;
    178,31,36;
    %    33, 144, 140; % microglia color    
    252,187,131]/255; 

figure;
for i=1:2
    currpcs = pcaH1_bc.Hlist{i};
    scatter(currpcs(1,:),currpcs(2,:),8,[0.5,0.5,0.5],'filled','MarkerFaceAlpha',0.5); hold on;
end

% now plot from E18 brain
i=3
for j=1:gmmdata.gmmlist{i}.NumComponents
    currpcs = pcaH1_bc.Hlist{i};
    currclusts = gmmdata.clusterlist{i};
    currcells = find(currclusts==j);
    currcol = collist(j,:);
    scatter(currpcs(1,currcells),currpcs(2,currcells),8,currcol,'filled','MarkerFaceAlpha',0.5); hold on;
end
axis tight
set(gca,'XTick',[0,5,10,15])
set(gca,'YTick',[0,5,10,15])
set(gca,'FontSize',16)
set(gcf,'pos',[360   519   310   262]);
print('figures/Fig1_E18-subpops.png','-dpng','-r200');


%% -------------------------------------------------------------
% Fig1d Make gene violin plots for markers
% -------------------------------------------------------------

%let's just compare Pax6 across E18 samples
df3 = full(log1.logdflist{3});

genes = {'PAX6','SOX2','NEUROG2','SOX9','ALDH1L1','GLUL','MAP2','TUBB3'};
[r,c] = findsubplotsize(length(genes));
figure;
for j=1:length(genes)
    gene = genes{j};
    gidx = findgeneinds(logdata.genes,{gene})
    idx1 = find(gmmdata.clusterlist{3}==2);
    idx2 = find(gmmdata.clusterlist{3}==5);
    idx3 = find(gmmdata.clusterlist{3}==1);
    idx4 = find(gmmdata.clusterlist{3}==3);

    dflist = {df3(gidx,idx1),df3(gidx,idx2),df3(gidx,idx3),df3(gidx,idx4)};
%    subplot(length(genes),1,j);
    subplot(r,c,j);
    violinPlot(dflist,'showMM',2,'divFactor',2);
    title(gene)
    ylim([-0,max(df3(gidx,:))])
end

set(gcf,'pos',[1881         448         654         233])
print('figures/Fig1_markers.pdf','-dpdf','-r300')
plotGeneViolin(gmmdata, pcaHdata_bc, logdata,'figures/genes/',brewermap(11,'Set2'), {'PAX6','TUBB3'}, {[3,2],[3,5],[3,1],[3,3]});%,[10,3]})

%% -----------------------------------------
% FIGURE 1g: PLOT HEATMAP FOR E18 CELLS
%------------------------------------------
ind = 3;
% load('/Users/sisichen/Dropbox (Thomson Lab)/Papers/RefMap/Analysis/baseline-pca/gmms/gmms_10_kmeans/6_SVZ-D0/SVZ-D0_8_gmm.mat')
% gmmdata.gmmlist{ind} = gmfitall;
% gmmdata.clusterlist{ind} = gmfitall.cluster(pcaH1.Hlist{ind}');

c = gmmdata.clusterlist{ind};
df = log1.logdflist{ind};
df = df(progidx,:);

% reorder the populations by sorting mu
mu = gmmdata.gmmlist{ind}.mu;
numcomp = size(mu,1);
[mu_r sortidx ] = sortrows(mu, 'd');
% ignore sorting
sortidx = 1:length(sortidx);

% remove microglia
%sortidx = sortidx([1,3:5])

comb_idx = [];
for i=1:length(sortidx)
    % instead of finding cells in pop#= i, use pop# = sortidx(i)
    currpop = sortidx(i);
    curridx = randsample(find(c==currpop),200,'true');
%     curridx = find(c==currpop)';
    comb_idx = [comb_idx,curridx];
end

df2 = scaleandsaturate(df(:,comb_idx),8,0);
figure; imagesc(df2)
g = caxis; 
set(gca,'YTick', 1:length(progidx), 'YTickLabel', proggenes);
set(gcf,'pos',[303   282   665   532]);
print(['figures/Fig1_',log1.dfnames{ind},'_prog_heatmap.pdf'], '-dpdf', '-r300')


%% -----------------------------------------
% FIGURE SI3: PLOT HEATMAP FOR SVZ cells
%------------------------------------------
ind = 4;
% load('/Users/sisichen/Dropbox (Thomson Lab)/Papers/RefMap/Analysis/baseline-pca/gmms/gmms_10_kmeans/6_SVZ-D0/SVZ-D0_8_gmm.mat')
% gmmdata.gmmlist{ind} = gmfitall;
% gmmdata.clusterlist{ind} = gmfitall.cluster(pcaH1.Hlist{ind}');

c = gmmdata.clusterlist{ind};
df = log1.logdflist{ind};
df = df(progidx,:);

% reorder the populations by sorting mu
mu = gmmdata.gmmlist{ind}.mu;
numcomp = size(mu,1);
[mu_r sortidx ] = sortrows(mu, 'd');

% remove microglia
%sortidx = sortidx([1,3:5])

comb_idx = [];
for i=1:length(sortidx)
    % instead of finding cells in pop#= i, use pop# = sortidx(i)
    currpop = sortidx(i);
    curridx = randsample(find(c==currpop),200,'true');
%     curridx = find(c==currpop)';
    comb_idx = [comb_idx,curridx];
end

df2 = scaleandsaturate(df(:,comb_idx),8,0);
figure; imagesc(df2)
g = caxis; 
set(gca,'YTick', 1:length(progidx), 'YTickLabel', proggenes);
set(gcf,'pos',[303   282   665   532]);
colorbar('southoutside')
print(['figures/SIFig3_',log1.dfnames{ind},'_prog_heatmap.pdf'], '-dpdf', '-r300')

%% -----------------------------------------
% Fig1f Celltype proportions
%------------------------------------------
% These numbers are assembled from data: 

prop = [71.4,13.8,12.4;67,16.2,16.8]/100;
figure;bar(prop',1,'LineStyle','none')
ylabel('proportions')
set(gca,'XTIckLabel',{'CN','NP','AP'})
set(gca, 'Fontsize',16);
legend('E18', 'SVZ')
set(gcf,'pos',[440   615   236   183])
print('figures/Fig1_proportions.pdf','-dpdf','-r300')


%% -----------------------------------------
% Build models for other datasets
%------------------------------------------
dfnames = { 'AdCo','P4', 'E18', 'SVZ-D0', 'SVZ-D3', 'SVZ-D5', 'SVZ-D7', 'SVZ-D18', ...
    'NSC-CTRL','NSC-CTRL-DIFF',...
    'SVZ-D5-CHIR', 'SVZ-D5-PDGF','SVZ-D5-BMP4',...
     'SVZ-D0-DIFF'};
rawdata = loadraw('../data/', dfnames, 'SVZ-D0','alldata','genes.tsv','mtx');


rawmult3 = loadscreen('../data/', 'NPC-MULT-3', 'Growth_BMP-', 'mult-3','genes.tsv','mtx');
rawmult4 = loadscreen('../data/', 'NPC-MULT-4', 'Diff_BMP-', 'mult-4','genes.tsv','mtx');
rawmult = mergestructures(rawmult3,rawmult4)

% add Growth_BMP- and Diff_BMP- to the data structure
rawdata.dflist = [rawdata.dflist,rawmult.dflist([1,3])];
rawdata.dfnames = [rawdata.dfnames, {'Growth_BMP-','Diff_BMP'}]
rawdata.sampnums = 1:17;

% load rawdata.mat
logdata = lognorm(rawdata);
genes = logdata.genes;

%% Remove RBCs: 
[logdata1, rawdata1] = removeRBCs(logdata, rawdata, 'mouse');

%% Filter data based on E18 P4, SVZ-D0 data. 
hvgthresh = 1.4;
filtdata = filterdata(logdata, rawdata, 'method', 'hvg', 'thresh', hvgthresh, 'removeRibo', true);


%% Project data using existing PC
pcaHdata = projectPCAdata(pcafeats, filtdata);

%% do the daisy chain CPD correct
% dflist = cellfun(@(x) x', pcaHdata.Hlist, 'Uniform', false);
% [newdflist, transforms] = daisy_CPD(dflist,500, 1,pcaHdata.dfnames);

% apply transforms to the growth and diff data: 
dflist2 = cellfun(@(x) x', pcaHdata.Hlist, 'Uniform', false);
 [newdflist2, transforms2] = daisy_CPD(dflist([8,16:17]),500, 1,pcaHdata.dfnames);

Rlist = [transforms.Rlist(1:8),transforms2.Rlist];
tlist = [transforms.tlist(1:8),transforms2.tlist];
 
newdflist3 = {}
dflist3 = dflist2([1:8,16:17])

for i=2:length(dflist3)
    currX = dflist3{i};
    [M,D2]=size(currX);
    
    % get list of transforms to apply:
    Torder = fliplr(1:(i-1));
    
    % Apply transformations backwards: 
    for j=1:length(Torder)
        Tnum = Torder(j);
        currR = Rlist{Tnum};
        currt = tlist{Tnum};
        newX = currX*currR' + ones(M,1)*currt';
        currX = newX;

    end
    newdflist3{i} = currX;
end

newdflist_final  = [newdflist(1:15),newdflist3{9:10}];
 

%% Make new pcaH1
pcaHdata_bc = pcaHdata;
pcaHdata_bc.Hlist = cellfun(@(x) x',newdflist_final,'Unif',false);

%% -------------------------------------------------------------
% SI fig? Scatter plot of all timepoints, corrected
% -------------------------------------------------------------
figure;
collist = flipud(viridis(15));
% combine all data and make a color vector
axlist=[];
for i=1:length(pcaHdata_bc.Hlist)
    axlist(i) = subplot(3,5,i);
    currpcs = pcaHdata_bc.Hlist{i};
    scatter(currpcs(1,:),currpcs(2,:),8,collist(i,:),'filled','MarkerFaceAlpha',0.5);hold on; 
    title(dfnames{i});
end
linkaxes(axlist,'xy')

%% Build gmms
rlist = 1:10;
gmmdata = buildgmm(pcaHdata_bc, 1000, 1:10, 1,4);

%% -------------------------------------------------------------
% Fig1h Plot model PCA renderings
% -------------------------------------------------------------
plotmodelPCA(gmmdata, pcaHdata_bc, 1:2, 100,'PCA', true, 'Fig1_model_renderings', true)
plotmodelPCA(gmmdata, pcaHdata, 1:2, 100,'PCA', true, 'Fig1_model_renderings_nonbatchcorrected', true)

%% -------------------------------------------------------------
% Fig1i-k proportions over time
% -------------------------------------------------------------

% make timecourse only data
pcaHdata_timecourse = pcaHdata_bc;
pcaHdata_timecourse.dfnames = pcaHdata_timecourse.dfnames(4:9)
pcaHdata_timecourse.Hlist = pcaHdata_timecourse.Hlist(4:9)
pcaHdata_timecourse.sampnums = pcaHdata_timecourse.sampnums(4:9)

gmmtimecourse = buildgmm(pcaHdata_timecourse, 1000, 1:5, 20,4);

% Check the heatmap for 
plotGMMheatmaplow(gmmtimecourse.clusterlist{13}, logdata.logdflist{4}, logdata.genes, true, progidx, false)  
% determine that 1,2,and 5 are committed neurons

alignedPops = popAlign(gmmtimecourse, pcafeats, 'bif', 13);
abmat = consolidate_w(alignedPops,gmmtimecourse);

% add rows 1,2,and 5 together
abmat2 = [abmat(3,:);abmat(4,:); sum(abmat([1,2,5],:))]

sampidx=gmmtimecourse.sampnums-3;
allmeans = [];
popnames = {'pre-astro','pre-neuro','comm. neuroblast'};
figure;
for refpop = 1:size(abmat2,1)
    subplot(size(abmat2,1),1,refpop)
    [means,stds]=errorplot(sampidx,abmat2(refpop,:),unique(sampidx),'LineWidth',3); hold on;
    allmeans = [allmeans;means];
    scatter(sampidx,abmat2(refpop,:),20,'k', 'filled', 'Jitter', 'on', 'JitterAmount', 0.1); hold on;
    set(gca,'XTIck', unique(sampidx),'XTickLabel', gmmtimecourse.dfnames, 'TickLabelInterpreter', 'none');
    set(gca,'FontSize',16)
    xtickangle(45)
    title(popnames{refpop});
    ylim([0,inf]);
    
end
set(gcf,'Position',[ 1930         -31         195         649])
print('figures/Fig1_abundances.pdf', '-dpdf', '-r300')

%%-------------------------------------
% Plot 2DProj of each dataset
%-------------------------------------
plot2dProj(gmmdata, pcaHdata_bc, 4, 1:10)

%%-------------------------------------
% Plot 2D errors for each datset
%-------------------------------------

[gerr,derr]= calcGMMerr(gmmdata, pcaHdata_bc, 1, 2000, binwidths, true);
print(['figures\SIFIg2a_Err_by_proj_index.pdf'],'-dpdf','-r300')

% Calculate binwidths for all dimensions
allH = [pcaHdata_bc.Hlist{:}];
ranges = range(allH');
binwidths = ranges./5;

 allgerrs = [];
 allderrs = [];
 stdgerr = [];
 stdderr = [];

 for i=1:length(gmmdata.gmmlist)

     currgerrs = [];
     currderrs = [];
     for j=1:10
         [gerr,derr]= calcGMMerr(gmmdata, pcaHdata_bc, i, 2000, binwidths, false);
         currgerrs = [currgerrs,gerr];
         currderrs = [currderrs,derr];
     end
     mgerr = mean(currgerrs);
     mderr = mean(currderrs);
     allgerrs = [allgerrs, mgerr];
     allderrs = [allderrs,mderr];
     
     stdgerr = [stdgerr, std(currgerrs)];
     stdderr = [stdderr, std(currderrs)];
 end

y = [allgerrs;allderrs]';
err = [stdgerr;stdderr]';
 
ngroups = size(y, 1);
nbars = size(y, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
 
figure; p=bar(y); hold on;
p(1).FaceColor=[0.8,0,0]
p(2).FaceColor='k';
ylim([0,20]);
for i=1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    e=errorbar(x,y(:,i),err(:,i)*0,err(:,i),'.', 'LineWidth',2); 
    e.Color='k';
end
ylabel('%error') 
set(gca,'FontSize',16)
set(gca,'Xtick',1:length(dfnames),'XTickLabel',dfnames)
xtickangle(45);
legend('model','data')
print(['figures\SIFig2_2D_error.pdf'],'-dpdf','-r300')

%% ===================================
%-------------------------------------
% Figure 2
%-------------------------------------

plotmodelPCA(gmmdata, pcaHdata_bc, 1:2, 100,'PCA', true, 'Fig2_model_renderings', true)

%=====================================
%%-------------------------------------
% check populations for SVZ-D0-DIFF
%-------------------------------------
plotGMMheatmaplow(gmmdata.clusterlist{15}, logdata.logdflist{15}, logdata.genes, true, [], true)  

%% -------------------------------------
% Plot heatmaps for all samples - unsupervised
%-------------------------------------
% 
% poplist = [2,4,9,10,15,16,17];
% for i=poplist
%     plotGMMheatmaplow([gmmdata.clusterlist{i}], [logdata.logdflist{i}], logdata.genes, true, progidx, false)
%     title(logdata.dfnames(i))
%     print(['Fig2_heatmaps_',logdata.dfnames{i},'.pdf'],'-dpdf','-r300');
% end

%% -------------------------------------
% Identifying genes from P4 data
% --------------------------------------
% 
% p4genes=plotGMMheatmaplow(gmmdata.clusterlist{2}, [logdata.logdflist{2}], logdata.genes, true,  [], true);
% p4idx = findgeneinds(logdata.genes, p4genes(1202:2373)); % will change if re-run
% plotGMMheatmaplow(gmmdata.clusterlist{2}, [logdata.logdflist{2}], logdata.genes, true,  p4idx, false)
% 
% p4idx2 = p4idx([298:360,927:1010,838:930,369:430,63:147,150:296]); % these are all the indexes for all cell types in the P4 brain
% p4idx3 = p4idx([298:360,927:1010,848:930]); %only neuronal, astrocytic and mitotic genes
% 
% % reorder p4idx2 genes: 
% normvals = max(logdata.logdflist{2}');
% p4_all_genes = plotGMMheatmaplow(gmmdata.clusterlist{2}, [logdata.logdflist{2}], logdata.genes, true,  p4idx2, true, normvals)
% p4idx_all = findgeneinds(logdata.genes,p4_all_genes);
% 
% p4idx_ANM = p4idx3; % only astrocyte neuron mitosois
% p4_ANM_genes = logdata.genes(p4idx2);
% 
% save p4genes.mat p4idx_all p4idx_ANM p4_all_genes p4_ANM_genes

%% -------------------------------------
% Fig 2c,f comparison heatmaps
% --------------------------------------

% re-cluster genes according to SVZ-D0-DIFF and NPC-D5-DIFF data 
compgenes = plotGMMheatmaplow([gmmdata.clusterlist{15};gmmdata.clusterlist{17}+max(gmmdata.clusterlist{15})], [logdata.logdflist{15},logdata.logdflist{17}], logdata.genes, true,  p4idx3, true,normvals)
set(gcf,'pos', [1760  149   678 428])

compgenes2 = compgenes([31:50,81:110,140:158]);
compinds2 = findgeneinds(logdata.genes,compgenes2);

normvals = max(logdata.logdflist{2}');

plotGMMheatmaplow(gmmdata.clusterlist{2}, [logdata.logdflist{2}], logdata.genes, true,  compinds2, false, normvals,[4])
set(gcf,'pos',[2745         347         111         420]);
print('Fig2_heatmaps_compgenes_P4_neuron.pdf','-dpdf','-r300')
plotGMMheatmaplow(gmmdata.clusterlist{2}, [logdata.logdflist{2}], logdata.genes, true,  compinds2, false, normvals,[8])
set(gcf,'pos',[2745         347         111         420]);
print('Fig2_heatmaps_compgenes_P4_astrocyte.pdf','-dpdf','-r300')
plotGMMheatmaplow(gmmdata.clusterlist{2}, [logdata.logdflist{2}], logdata.genes, true,  compinds2, false, normvals,[10])
set(gcf,'pos',[2745         347         111         420]);
print('Fig2_heatmaps_compgenes_P4_AP.pdf','-dpdf','-r300')
plotGMMheatmaplow(gmmdata.clusterlist{15}, [logdata.logdflist{15}], logdata.genes, true,  compinds2, false,normvals)
set(gcf,'pos', [ 2710        -147         347         420])
print('Fig2_heatmaps_compgenes_SVZ-D0-DIFF.pdf','-dpdf','-r300')
plotGMMheatmaplow(gmmdata.clusterlist{17}, [logdata.logdflist{17}], logdata.genes, true,  compinds2, false,normvals,[1,3])
set(gcf,'pos', [ 2538         347         191         420])
print('Fig2_heatmaps_compgenes_NPC-D5-DIFF.pdf','-dpdf','-r300')

plotGMMheatmaplow(gmmdata.clusterlist{2}, [logdata.logdflist{2}], logdata.genes, true,  compinds2, false, normvals)
set(gcf,'pos',[2745         347         111         420]);
title('P4 all')
print('Fig2_heatmaps_compgenes_P4_all.pdf','-dpdf','-r300')
plotGMMheatmaplow(gmmdata.clusterlist{15}, [logdata.logdflist{15}], logdata.genes, true,  compinds2, false, normvals)
set(gcf,'pos',[2745         347         111         420]);
title('SVZ-D0-DIFF all')
print('Fig2_heatmaps_compgenes_SVZ-D0-DIFF_all.pdf','-dpdf','-r300')

plotGMMheatmaplow(gmmdata.clusterlist{16}, [logdata.logdflist{16}], logdata.genes, true,  compinds2, false, normvals)
set(gcf,'pos',[2745         347         111         420]);
title('NPC-D5-GROWTH all')
print('Fig2_heatmaps_compgenes_NPC-D5-GROWTH_all.pdf','-dpdf','-r300')

plotGMMheatmaplow(gmmdata.clusterlist{17}, [logdata.logdflist{17}], logdata.genes, true,  compinds2, false, normvals)
set(gcf,'pos',[2745         347         111         420]);
title('NPC-D5-DIFF all')
print('Fig2_heatmaps_compgenes_NPC-D5-DIFF_all.pdf','-dpdf','-r300')


% %-------------------------------------
% % Run PopAlign (Can use this information later to generate dendrogram)
% %-------------------------------------

alignedPops = popAlign(gmmdata, pcafeats, 'test2ref', 2)
A=alignedPops.alignidx

% Find the alignments specifically for 
diffinds = [find(A.test_idx ==17);find(A.test_idx ==15)]
JDlist = A.JD(diffinds)


alignedPops2 = popAlign(gmmdata, pcafeats, 'ref2test', 2)
A2 = alignedPops2.alignidx

%%
%-------------------------------------
% Fig 2 Extra Diff scatterplot of SVZ-D0-DIFF
%-------------------------------------

bkgd = [pcaHdata_bc.Hlist{1:3}];
bkgdidx = randsample(size(bkgd,2),10000,1);
sampbkgd = bkgd(:,bkgdidx);
poplist = [4,9,10,15,16,17];

for i=[7,8,9,10,15,16,17]% i=[9,10,15,16,17];

    figure; 
    scatter(sampbkgd(1,:),sampbkgd(2,:),10,[0.5,0.5,0.5],'filled','MarkerFaceAlpha',0.5);hold on; 
    % now plot from SVZ-D0-DIFF 
    collist = brewermap(gmmdata.gmmlist{i}.NumComponents ,'Set2')

    for j=1:gmmdata.gmmlist{i}.NumComponents
        currpcs = pcaHdata_bc.Hlist{i};
        currclusts = gmmdata.clusterlist{i};
        currcells = find(currclusts==j);
        currcol = collist(j,:);
        scatter(currpcs(1,currcells),currpcs(2,currcells),8,currcol,'filled','MarkerFaceAlpha',0.5); hold on;
    end
    %title(gmmdata.dfnames{i}, 'Interpreter','none')
    axis tight
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    set(gca,'FontSize',16)
    set(gcf,'pos',[360   519   310   262]);
     %print(['Fig2_scatterplots_',logdata.dfnames{i},'.png'],'-dpng','-r300');    

end

%-
