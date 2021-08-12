function [newlogdata, newrawdata] = removeRBCs(logdata, rawdata, species)
% function [newlogdata, newrawdata] = removeRBCs(logdata, rawdata, species)
% This function removes RBCs from log normalized data based on RBC marker
% genes. The thresholds are calculated on the log normalized data but both
% logdata and rawdata structures are filtered. 
 
nfiles = length(logdata.logdflist);

mouseRBCgenes = {'HBB-BT','HBB-BS' ,'HBA-A1','HBA-A2'};
humanRBCgenes = {'HBB', 'HBA1', 'HBA2'};

switch species
    case 'mouse'
        genelist = cellfun(@(x) find(strcmp(logdata.genes, x)), mouseRBCgenes);
    case 'human'
        genelist = cellfun(@(x) find(strcmp(logdata.genes, x)), humanRBCgenes);
end

% Sum the RBC genes for each cell and keep track of the sum
allRBCsum=[];
RBCsumlist = cell(1,nfiles);
for i=1:nfiles
    currdata = logdata.logdflist{i};
    RBCsum = sum(currdata(genelist,:));
    allRBCsum = [allRBCsum,RBCsum];
    RBCsumlist{i} = RBCsum;
end

% find the threshold;
T = multithresh(full(allRBCsum));

% Loop over datasets and remove RBCs from each log normalized dataset
newlogdflist = cell(1,nfiles);
newrawdflist = cell(1,nfiles);

for i=1:nfiles
    currRBCsum = RBCsumlist{i};
    keepinds = find(currRBCsum<T);
    
    % filter log data
    currlogdata = logdata.logdflist{i};
    newlogdflist{i} = currlogdata(:,keepinds);
    
    % filter raw data
    currrawdata = rawdata.dflist{i};
    newrawdflist{i} = currrawdata(:,keepinds);

end

newlogdata = logdata;
newlogdata.logdflist = newlogdflist;

newrawdata = rawdata;
newrawdata.dflist = newrawdflist;
