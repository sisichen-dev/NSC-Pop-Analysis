function Hdata = projectdata(feats, filtdata)
% function Hdata = projectdata(feats, filtdata);
% returns structure with filtered data
% Hdata.expname
% Hdata.dnames
% Hdata.refname
% Hdata.Hlist = {H1, H2, ...}
% Hdata.W = featdata.W;

fprintf('\n')
disp('************************************')
disp('Projecting data into feature space....')
disp('************************************')
fprintf('\n')

% Gathering 
filtdatalist = filtdata.filtdflist;
dfnames = filtdata.dfnames;
Hfolder = 'Hdata/';
% Hnames = strcat(Hfolder, cellstr(num2str(filtdata.sampnums'))','_', dfnames, '_',num2str(feats.k),'_Hdata.mat');
Hnames = strcat(Hfolder, dfnames, '_Hdata.mat');

expname = filtdata.expname;
W = feats.W;

% accumulate indexes of files that have already been filtered
rem_inds = [];

% Check if Hfolder already exists, if not make it.
if (exist(Hfolder,'dir')==7)
    % Remove files which already exist in filt folder
    for i=1:length(filtdatalist)
        if (exist(Hnames{i},'file'))==2
            disp([Hnames{i}, ' already exists']);
            rem_inds = [rem_inds,i];
        end
    end
else
    mkdir(pwd,Hfolder)
end

proj_inds = setdiff(1:length(Hnames),rem_inds);

% Calculate H for all datasets
for i = proj_inds
    currdata = filtdatalist{i};
    disp(['Projecting dataset into H space : ', Hnames{i}]);
    currH = nnls(W,currdata');
    fprintf('\n');
    Hlist{i}=currH;
    % save H in its own file
    save(Hnames{i}, 'currH');
end

% Load in H for all datasets
Hlist = {};
for i=1:length(Hnames)
    
    disp(['Loading in: ',Hnames{i}])
    load(Hnames{i});
    Hlist{i}=currH;
end


% separate the test and reference files
n = length(filtdata.refname);

Hdata.expname = filtdata.expname;
Hdata.dfnames = filtdata.dfnames;
Hdata.refname = filtdata.refname;
Hdata.sampnums = filtdata.sampnums;
Hdata.Hlist = Hlist;
Hdata.W = feats.W;