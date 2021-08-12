function plotGMMheatmap(gmmdata, logdata, ind)
% function plotGMMheatmap(gmmdata, logdata, ind)
%  This function takes a sample index and then plots a gene expression
%  heatmap 

currgmm = gmmdata.gmmlist{ind};
genes = logdata.genes;
currcluster = gmmdata.clusterlist{ind};
currlogdata = logdata.logdflist{ind};

plotGMMheatmaplow(currgmm, currcluster, currlogdata, genes)

