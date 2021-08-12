function rawdata = cullsamples(rawdata,minN)
% function rawdata = cullsamples(rawdata,minN)
% remove all samples which have fewer than minN number of cells

cellnums=cellfun(@(x) size(x,2), rawdata.dflist);
remIdx = find(cellnums<minN);
newmeta = rawdata.meta;
% find all indices of the samples to be removed and remove them from the
% metadata
for i=1:length(remIdx)
    idx = remIdx(i);
    meta_remIdx = find(newmeta{:,2}==idx);
    newmeta(meta_remIdx,:)=[];
end

% remove these samples from the dflist
rawdata.dflist(remIdx) = [];
rawdata.dfnames(remIdx) = [];
rawdata.meta = newmeta;
rawdata.sampnums(remIdx) = [];
