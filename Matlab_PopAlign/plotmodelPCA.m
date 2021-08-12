function plotmodelPCA(gmmdata, Hdata, dims, npixel,feattype, weighted, fname, plotTogether, cmap)
% plotmodelPCA(gmmdata, Hdata, dims, npixel,feattype, weighted, fname, plotTogether)
% plots PCA renderings of models
% dims is a vector indicating which dimensions will be plotted

figfolder = 'figures/';
orange = [1.0000 , 0.4235 , 0.0471];
gray = [0.3, 0.3, 0.3];

if nargin==3
    npixel=100;
    feattype = 'PCA';
    plotTogether=true;
    cmap = 'jet';
elseif nargin==4
    plotTogether=true;
    cmap = 'jet';
elseif nargin ==5
    cmap = 'jet';
end

% Check if gmmfolder already exists, if not make it.
if (exist(figfolder,'dir')==7)
    
else
    mkdir(pwd,figfolder)
end

% make name depending on whether weighted or not: 
if weighted
    fname2 = [fname, '_weighted']
else
    fname2 = [fname, '_unweighted']
end

% Initializing variables
% keep only the best models
modellist = gmmdata.gmmlist(gmmdata.bestlist);
totaln = length(modellist)
nameslist = gmmdata.modelnames(gmmdata.bestlist);

Hlist = Hdata.Hlist;

%% Plot model rendering in PC space
fprintf('\n')
disp('************************************')
disp('Plotting models in PCA space')
disp('************************************')
fprintf('\n')

alldata = [Hdata.Hlist{1:end}]';
N = min(5000,size(alldata,1));

if strcmp(feattype,'ONMF')
    sampdata = full(alldata(randsample(size(alldata,1),N),:));
    [coeff,score, LATENT, TSQUARED, EXPLAINED, pcamean] = pca(sampdata);
elseif strcmp(feattype,'PCA')
    coeff = eye(size(alldata,2));
    pcamean = zeros(1,size(alldata,2));
end
allscores = (alldata-pcamean)*coeff;

% define space 
alpha = 0.15; % size of margin
xmin = min(allscores(:,dims(1)))-alpha;
xmax = max(allscores(:,dims(1)))+alpha;
ymin = min(allscores(:,dims(2)))-alpha;
ymax = max(allscores(:,dims(2)))+alpha;

delx = (xmax - xmin)/(npixel-1);
dely = (ymax - ymin)/(npixel-1);
X = xmin:delx:xmax;
Y = ymin:dely:ymax;

% %%%%%%%%%%%%%%
% % solve linear system of equations to find transformation between PCA and
% % pixels
% 
% syms b0 b1 a11 a12 a21 a22
% 
% eqn1 = -b0 == a11* xmin + a12 * ymin;
% eqn2 = npixel - b1 == a21 * xmin + a22 * ymin;
% eqn3 = -b0 == a11*xmin + a12 * ymax;
% eqn4 = -b1 == a21* xmin + a22 * ymax;
% eqn5 = npixel - b0 == a11 * xmax + a12 * ymax;
% eqn6 = -b1 == a21 * xmax + a22 * ymax;
% eqn7 = npixel - b0 == a11 * xmax + a12 * ymin;
% eqn8 = npixel - b1 == a21 * xmax + a22 * ymin;
% 
% [A,B] = equationsToMatrix([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8], [b0,b1,a11,a12,a21,a22]);
% 
% Z = linsolve(A,B);
% Zreal = eval(Z);
% 
% b = [Zreal(1);Zreal(2)];
% A = [Zreal(3),Zreal(4);Zreal(5),Zreal(6)];

%%%%%%%%%%%%%%
figure;
[r,c]=findsubplotsize(totaln);

for j=1:(totaln)
    % project model down into PCA space
    currmodel = modellist{j};
    currname = nameslist{j};
    modelmus = currmodel.mu;
    muproj = (modelmus - pcamean) * coeff ;
    newmus = muproj(:,dims);
    nmu = size(muproj,1);

    newSigma=ones(2,2,nmu);
    for i=1:nmu
        currSigma = currmodel.Sigma(:,:,i);
        projSigma = coeff'* currSigma *coeff;
        newSigma(:,:,i) = projSigma(dims,dims);
    end
    
    % compute likelihoods for each model separately: 
    allprobs=zeros(length(X),length(Y));
    for g = 1:nmu
        subgm = gmdistribution(newmus(g,:),newSigma(:,:,g));
        
        probs=zeros(length(X),length(Y));
        for i=1:length(X)
            x1=X(i);
            for k=1:length(Y)
                y1=Y(k);
                prob = subgm.pdf([x1,y1]);
                probs(k,i) = prob;
            end
        end
        probs = probs ./ (sum(sum(probs)));
        if weighted
            allprobs = allprobs + probs*currmodel.ComponentProportion(g);
        else % not weighted
         allprobs = allprobs + probs;
        end
    end
    
    mu_im_x = (newmus(:,1)-xmin)/(xmax-xmin)*(npixel);
    mu_im_y = npixel-(newmus(:,2)-ymin)/(ymax-ymin)*(npixel);
    
    pdel = npixel/10;
    
    if plotTogether
        subplot(r,c,j) 
        multfactor = 100;
    else
        multfactor = 1000;
        figure;
    end

    colormap(jet);
    imagesc(flipud(log(allprobs*10000+1))); shading('interp'); hold on;
    scatter(mu_im_x, mu_im_y, currmodel.ComponentProportion*multfactor, 'k', 'filled', 'MarkerFaceAlpha',0.5)
    text((mu_im_x+5), mu_im_y + 0.02, strsplit(num2str(1:nmu)), 'Color','k', 'FontSize', 18)
    
    % if want to reverse x
%     imagesc(fliplr(flipud(log(allprobs*10000+1)))); shading('interp'); hold on;
%     scatter(npixel-mu_im_x, mu_im_y, currmodel.ComponentProportion*multfactor, 'k', 'filled', 'MarkerFaceAlpha',0.5)
%     text(npixel-(mu_im_x+5), mu_im_y + 0.02, strsplit(num2str(1:nmu)), 'Color','k', 'FontSize', 18)    
     
    title(currname, 'interpreter','none', 'FontSize', 12)
    %set(gca,'XTick', 1:pdel:npixel, 'XTickLabel', round(X(1:pdel:npixel),2))
    %set(gca,'YTick', 1:pdel:npixel, 'YTickLabel', fliplr(round(Y(1:pdel:npixel),2)))
    set(gca,'YTick', []);
    set(gca, 'XTick', []);
   
    %xlabel(['PC',num2str(dims(1))], 'FontSize',12)
    %ylabel(['PC',num2str(dims(2))], 'FontSize',12)
    
    if ~plotTogether
        set(gcf,'Position',[0,0,200,150])
        print([figfolder, fname2,'_',num2str(j),'_',currname,'.pdf'],'-dpdf','-r300')        
    end
    
    disp(['plotted pca figure for: ', num2str(i), ' ', currname,]);
   
end

if plotTogether
    set(gcf,'Position',[0,0,200*c,150*r])
    print([figfolder, fname2,'.pdf'],'-dpdf','-r300') 
end