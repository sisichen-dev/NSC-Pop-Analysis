function [pcafeats, pcaHdata] = makePCAfeat(filtdata, sampsize,ndim)
% function [pcafeats, pcaHdata] = makePCAfeat(filtdata, sampsize,ndim)
%
% Uses PCA to find featureset for all 
% ndim is the number of dimensions used for calculations

dfnames = filtdata.dfnames;
nfiles = length(dfnames);
alldata = [filtdata.filtdflist{:}];
N = min(sampsize,size(alldata,2));
sampdata = full(alldata(:,randsample(size(alldata,2),N)))';
[V, T, LATENT, TSQUARED, EXPLAINED, pcamean] = pca(sampdata);

%% Now store all PCA values in a structure
pcafeats.W = V(:,1:ndim);
pcafeats.T = T;
pcafeats.pcamean = pcamean;
pcafeats.latent = LATENT;
pcafeats.explained = EXPLAINED;
pcafeats.feattype = 'PCA';
pcafeats.k = ndim;

%% also store all the projected data into a new Hdata structure: 

Hlist={};
for i=1:nfiles
    temp = (filtdata.filtdflist{i} - pcamean')' * V;
    % store full pca scores in pcafeats: 
    pcafeats.pclist{i} = temp;
    Hlist{i}=temp(:,1:ndim)';
end

pcaHdata.expname = filtdata.expname;
pcaHdata.dfnames = filtdata.dfnames;
pcaHdata.refname = filtdata.refname;
pcaHdata.Hlist = Hlist;
pcaHdata.W = V;
pcaHdata.sampnums = filtdata.sampnums;

