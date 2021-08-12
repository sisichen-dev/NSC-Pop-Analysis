function logdata = lognorm(rawdata, alpha)
% logdata = lognorm(rawdata, alpha)
% returns structure with log-normalized data
% logdata.expname = experiment name;
% logdata.dfnames = {dfname1, dfname2, ...}
% logdata.logdflist = {logdf1, logdf2, ...}
% logdata.refname = reference name
% logdata.genes = genes;
% logdata.sampnums = sample numbers | for screens number correlates to
% sampnum in the metadata file

if nargin==1
    alpha = 3000 ; %default value for alpha
end

dfnames = rawdata.dfnames;
nfiles=length(dfnames);
logfolder = 'logdata/';
% lognames = strcat(logfolder,cellstr(num2str(rawdata.sampnums'))','_',dfnames, '_logdata.mat'); % to log
lognames = strcat(logfolder,dfnames, '_logdata.mat'); % to log


fprintf('\n')
disp('************************************')
disp('Normalizing and logging data....')
disp('************************************')
fprintf('\n')

% accumulate indexes of files that have already been normalized
rem_inds = [];

% Check if logfolder already exists, if not make it.
if (exist(logfolder,'dir')==7)
    
    % Remove files which already exist in log folder
    for i=1:nfiles
        if (exist(lognames{i},'file'))==2
            disp([lognames{i}, ' already normalized.'])
            rem_inds = [rem_inds, i];
        end
    end
else
    mkdir(pwd,logfolder)
end

% get indices of files that need to be log-normalized
tolog_inds = setdiff(1:length(lognames), rem_inds);

% Iterate through indices and log normalize data 
for i = tolog_inds
    
    disp(['Normalizing: ',lognames{i}])
    
    % Load file depending on type
    currdf = rawdata.dflist{i};
    
    % normalize each column by total number of transcripts
    tot = sum(currdf);
    currdf2 = currdf./tot * alpha;

    % take the log + 1;
    logdf = log1p(currdf2);
   
    % save data
    save(lognames{i}, 'logdf')
   
end

% load in logdflist from all files
logdflist={};

for i=1:length(lognames)
    
    disp(['Loading in: ',lognames{i}])
    load(lognames{i});
    logdflist{i}=logdf;
end

% Store the logged data

logdata.expname = rawdata.expname;
logdata.dfnames = dfnames;
logdata.logdflist = logdflist;
logdata.refname = rawdata.refname;
logdata.genes = rawdata.genes;
logdata.sampnums = rawdata.sampnums;

