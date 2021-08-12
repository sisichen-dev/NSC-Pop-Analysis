function rawdata = loadscreen(foldername, dfname, refname, expname, gname, ftype)
% function rawdata = loadscreen(foldername, dfname, refname, expname, gname, ftype)
% returns structure with raw data
% rawdata.dflist = {df1, df2, ...}
% rawdata.dfnames = {name1, name2 , ...}
% rawdata.genes = genes
% rawdata.expname = expname;
% rawdata.refname = reference sample name
% rawdata.sampnums = 1:length(dfnames)
% rawdata.meta = metadata file

% foldername is normally: 'data/'

% Read in gene names
test = tdfread([foldername,gname]);
gene1 = fieldnames(test);
gene1 = gene1{2};
mostgenes = upper(cellstr(strtrim(getfield(test,gene1))));
rawdata.genes = {gene1, mostgenes{:}}';

fprintf('\n')
disp('************************************')
disp('Loading raw data from a multiplexed screen')
disp('************************************')
fprintf('\n')

% Read in metadata 
metafile = [foldername, dfname, '/', dfname, '-meta.csv'];
currmeta = readtable(metafile);

% Read in data 
currfile = [foldername, dfname, '/', dfname, '.',ftype];
if ftype=='mtx'
    df = mmread(currfile);
elseif ftype=='mat'
    load(currfile);
    df=data;
end

% segment data into samples based on metadata: 
samples = unique(currmeta.sample_number);
samples(isnan(samples)) = [];
numsamples = length(samples);
dflist = {};
dfnames = {};
varnames = currmeta.Properties.VariableNames(3:end);

for i=1:numsamples
    sampnum = samples(i);
    % extract the data from those samples: 
    curr_idx = find(currmeta.sample_number==sampnum);
%   
%     % for certain cases: 
%     currsample = samples{i};
%     curr_idx = find(strcmp(currmeta{:,2},currsample));
    
    % extract dataframe for the current index
    currdf = df(:,curr_idx);
    dflist = [dflist,{currdf}];   
 
    % if sample_id is already on the list, use it directly
    if (sum(contains(currmeta.Properties.VariableNames,'sample_id'))>0)
        name = currmeta.sample_id(curr_idx(1));
        
    elseif (size(currmeta,2)>2)       
        % Assemble a name based on the metadata entries
        curr_var_values = table2array(currmeta(curr_idx(1),3:end));
        if isnumeric(curr_var_values)
            keepvars = find(curr_var_values>0);

            if isempty(keepvars)
                name = 'control';
            else
                curr_var_values_keep = strsplit(num2str(curr_var_values(keepvars)));
                name = strcat(varnames(keepvars), '-', curr_var_values_keep);
                name = strjoin(name,'_');
            end
        else
            name = strjoin(curr_var_values,'_');
        end 
    
    else
        name = num2str(sampnum); % use sample number as name
    end

    dfnames = [dfnames,name];
end

% keep only the meta data that is assigned a sample number
% meta2 = currmeta((~isnan(currmeta.sample_number)),:);
meta2 = currmeta(currmeta.sample_number>0,:);

rawdata.dflist = dflist;
rawdata.dfnames = dfnames;
rawdata.refname = refname;
rawdata.expname = expname;
rawdata.meta = meta2;
rawdata.sampnums = 1:numsamples;


