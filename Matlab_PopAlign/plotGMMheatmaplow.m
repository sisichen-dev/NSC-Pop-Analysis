
function keepgenes = plotGMMheatmaplow(currcluster, currlogdata, genes, dispgenes, geneinds, clustgenes, normvals, clusterstoshow)
% function plotGMMheatmaplow(currcluster, currlogdata, genes, dispgenes, geneinds, clustgenes, normval, clusterstoshow)

%  This is a low level version of plotGMMheatmap. It takes the gmm, cluster
%  assignments list, original log data, and genes directly. 
% currcluster = list of cluster assignments for cells 
% currlogdata = logged data frame
% genes = gene names for currlogdata (must be the same length)
% dispgenes = boolean specificying whether gene names will be displayed
% geneinds = specific gene indices to use, can be [] if none specified
% clustgenes = boolean speficying whether genes should be clustered

if nargin==3
    dispgenes = true;
    geneinds = [];
elseif nargin==4
    geneinds = [];
end

if (length(genes)~=size(currlogdata,1))
    error('variable genes not the same length as the logdata')
end

if ~exist('clusterstoshow')
    clusterstoshow = unique(currcluster)';
end

clusterstoshow

numcells =[];
for i = clusterstoshow
    i
    currnum = sum(currcluster==i);
    numcells = [numcells, currnum];
end

sampnumber = min([100, numcells]);
%sampnumber = 100;

% make a sampled dataset: 
allkeep = [];
for i = clusterstoshow
    currcells = find(currcluster == i);
    currkeep = currcells(randsample(length(currcells), sampnumber)); 
    allkeep = [allkeep; currkeep];
end

clusterskeep = currcluster(allkeep);
[clustersorted, I] = sortrows(clusterskeep);
midtick = cellfun(@(x) mean(find(clustersorted==x)), num2cell(clusterstoshow));

sampdf = currlogdata(:,allkeep);

backsampdf = exp(sampdf')-1;
size(backsampdf)

% Now if no gene list is supplied, let's extract genes using the hvg
% threshold
if isempty(geneinds)
    geneinds = extractgenes(backsampdf', 'hvg',1.1,'figures/','NSC2', false)
 
%     % add all high-expression genes: 
%     highgeneinds = find(log10(mean(sampdf'))>(-0.5));
%     geneinds = unique([geneinds', highgeneinds]);

end

if isempty(normvals)
    normvals = max(sampdf');
end

sampdf2 = (sampdf' ./ normvals)';

% make the sampled heatmap: 
sampdf3 = sampdf2(geneinds,:);

if clustgenes
    [rowidx, colidx] = clusterMat(sampdf3, 'cor', 'single');
else
    rowidx = 1:length(geneinds);
end

%saturate the data
sampdf3(sampdf3>0.8) = 0.8;

% only cluster the genes: 
figure; 
ax1 = subplot(10,1,1:9)
imagesc(sampdf3(rowidx,:))
if dispgenes
    set(gca,'YTick', 1:length(geneinds), 'YTickLabel', genes(geneinds(rowidx)))
end
colorbar;
set(gca,'XTick',[]);
ax2=subplot(10,1,10);
imagesc(clustersorted(:)');
colormap(ax2,brewermap(length(unique(currcluster)),'Set2'))
set(gca,'YTick', [],'XTick', midtick, 'XTickLabel', num2cell(clusterstoshow))
colorbar

keepgenes = genes(geneinds(rowidx));