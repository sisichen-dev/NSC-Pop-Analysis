function filtdata = filterdata(logdata, rawdata, varargin)
% function filtdata = filterdata(logdata, rawdata, 'PARAM1', val1, 'PARAM2', val2, 'PARAM3', val3)
%
% This function filters the genes in the datasets
% 'method': must be 'mean' or 'hvg'
% 'thresh': 
%   for 'mean' thresholding: is the value of the gene mean
%   for 'hvg' thresholding:  log2(thresh) is the value of the offset from
%   the diagonal line for thresholding genes
% 'removeRibosomal': 
%   true: remove all genes starting with RPS or RPL
%   false: leaves genes as is
% returns structure with filtered data
% filtdata.reflist = {refdf1, refdf2, ...}
% filtdata.testlist = {testdf, testdf2, ...}
% filtdata.genelist = genelist

p = inputParser;

addRequired(p, 'logdata')
defaultmethod = 'mean';
checkString = @(s) any(strcmp(s,{'mean','hvg'}));
addParameter(p, 'method',defaultmethod,checkString);
addParameter(p, 'thresh', 0.01, @isnumeric);
addParameter(p, 'removeRibo', false, @islogical);

parse(p,'logdata',varargin{:});
method = p.Results.method;
thresh = p.Results.thresh;
removeRibo = p.Results.removeRibo;

fprintf('\n')
disp('************************************')
disp('Filtering data....')
disp('************************************')
fprintf('\n')

dfnames = logdata.dfnames;
nfiles = length(dfnames);
filtfolder = 'filtdata/';
% filtnames = strcat(filtfolder, cellstr(num2str(rawdata.sampnums'))','_', dfnames, '_filtdata.mat')
filtnames = strcat(filtfolder, dfnames, '_filtdata.mat')


% Check if filtfolder already exists, if not make it.
if ~(exist(filtfolder,'dir')==7)
    mkdir(pwd,filtfolder)
end

% Check to see if all datasets are in filter folder, if so, load all
% filtered dfs into a list
currfiles = dir([filtfolder, '*filtdata.mat']);
currnames = strcat(filtfolder,{currfiles.name});
if ~isequal(filtnames,currnames)
     filtdflist = filtall(logdata,rawdata, method, thresh, removeRibo);
else
    disp('All data already filtered.')
    % load in all filtered data iteratively
    filtdflist={};
    for i=1:length(filtnames)
        disp(['Loading in dataset: ', filtnames{i}])
        load(filtnames{i});
        filtdflist{i} = filtdf;
    end
end

load([filtfolder, 'allinds.mat'])

% Add fields to structure
filtdata.expname = logdata.expname;
filtdata.dfnames = logdata.dfnames;
filtdata.refname = logdata.refname;
filtdata.filtdflist = filtdflist;
filtdata.genes = logdata.genes(allinds);
filtdata.geneinds = allinds;
filtdata.sampnums = logdata.sampnums;

function filtdflist = filtall(logdata, rawdata,method, thresh, removeRibo)

dfnames = logdata.dfnames;
nfiles = length(dfnames);
filtfolder = 'filtdata/';
mkdir('figures', 'preprocessing');
filtfigfolder = 'figures/preprocessing/';
% filtnames = strcat(filtfolder,cellstr(num2str(rawdata.sampnums'))','_',  dfnames, '_filtdata.mat');
filtnames = strcat(filtfolder, dfnames, '_filtdata.mat');

logdflist = logdata.logdflist;
genes = logdata.genes;

% load in data from logfolder and get unique genes from each
switch(method)
    case 'mean'
        % Threshold the genes for each sample according to mean threshold
        dflist = logdata.logdflist;
    case 'hvg'
        dflist = rawdata.dflist;
end
alldf = [];
geneindlist={};

if (exist([filtfolder,'allinds.mat'])==2)
    disp('Loading existing gene list')
    load([filtfolder,'allinds.mat'])
else 
    disp('Extracting genelist from datasets' )
    
    % assemble combodf while looping:
    n = min(cellfun(@(x) size(x,2), dflist));
    s = min(n,1000);
    s = min(round(8000/nfiles), s); % take at maximum 8000 cells total
    
    for i=1:nfiles
        df=dflist{i}; % if method == hvg, this df will be raw
        currgeneinds = extractgenes(df, method, thresh, filtfigfolder, dfnames{i}, true);
        geneindlist{i} = currgeneinds;
        combodf{i} = df(:,randsample(size(df,2),s));
    end
    
    % find highly variable genes for the sampled combined dataframe and
    % also include these genes in the analysis
    combodf = [combodf{:}];
    currgeneinds = extractgenes(combodf, method, thresh, filtfigfolder, 'all', true);
    geneindlist{i+1} = currgeneinds;
    
    % Take the union of all
    allinds=geneindlist{1};
    for i=2:length(geneindlist)
        allinds = union(allinds, geneindlist{i});
    end
    
    % if remove ribosomal is true, remove all genes which begin with RPS
    % and RPL
    if removeRibo
        currgenes = genes(allinds);
        removeinds = find(startsWith(currgenes, {'RPL','RPS'}));
        allinds(removeinds) =[];
    end
    
end

disp(['Total Number of Genes: ', num2str(length(allinds))]);

% filter and save the data
filtdflist={};
for i=1:length(logdflist)
    disp(['Filtering dataset: ', dfnames{i}])
    
    currlogdf = logdflist{i};
    filtdf = currlogdf(allinds,:);
    filtdflist{i} = filtdf;
    save(filtnames{i}, 'filtdf')
end


save([filtfolder, 'allinds.mat'], 'allinds');

        
