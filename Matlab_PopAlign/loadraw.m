function rawdata = loadraw(foldername, dfnames, refname, expname, gname, ftype)
% function rawdata = loadraw(foldername, dfnames, refname, expname, gname, ftype)
% returns structure with raw data
% rawdata.dflist = {df1, df2, ...}
% rawdata.dfnames = {name1, name2 , ...}
% rawdata.genes = genes
% rawdata.expname = expname;
% rawdata.refname = reference sample name
% rawdata.sampnums = 1:length(dfnames);
% foldername is normally: 'data/'

nfiles=length(dfnames);

% Read in gene names
test = tdfread([foldername,gname]);
gene1 = fieldnames(test);
gene1 = gene1{2};
mostgenes = upper(cellstr(strtrim(getfield(test,gene1))));
rawdata.genes = {gene1, mostgenes{:}}';


fprintf('\n')
disp('***********************')
disp('Loading raw data ')
disp('***********************')
fprintf('\n')

rawdata.dflist = {};
% Iterate through files and log normalize data 
for i=1:length(dfnames)
    
    disp(['Loading in: ',dfnames{i}])
    
    % Load file depending on type
    currfile = [foldername, dfnames{i}, '.',ftype];
    
    if ftype=='mtx'
        currdf = mmread(currfile);
    elseif ftype=='mat'
        load(currfile);
%         currdf=array;
        currdf=data;
    end
    rawdata.dflist{i}=currdf;
end

rawdata.dfnames = dfnames;
rawdata.refname = refname;
rawdata.expname = expname;
rawdata.sampnums = 1:length(dfnames);



