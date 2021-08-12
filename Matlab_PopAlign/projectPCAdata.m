function pcaHdata = projectPCAdata(pcafeats, filtdata)
% function pcaHdata = projectPCAdata(pcafeats, filtdata)
% This function projects data into a pca space defined by the input
% pcafeats dataset 

numfiles = length(filtdata.filtdflist);

Hlist = {};

for j=1:numfiles

    currdf = filtdata.filtdflist{j};
    % First subtract the mean
    m1 = currdf - pcafeats.pcamean';

    % Then calculate projected data and rotate to match
    projdata = m1' * pcafeats.W;
    projdata = projdata';
    
    % Add projected data into Hlist
    Hlist{j} = projdata;
    
end

pcaHdata.expname = filtdata.expname;
pcaHdata.dfnames = filtdata.dfnames;
pcaHdata.refname = filtdata.refname;
pcaHdata.sampnums = filtdata.sampnums;
pcaHdata.Hlist = Hlist;
pcaHdata.W = pcafeats.W;
