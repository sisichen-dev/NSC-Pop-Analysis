function plot2dProj(gmmdata, Hdata, sampleInd, fplotinds)

gmmlist = gmmdata.gmmlist;%[{gmmdata.refmodel}, gmmdata.testlist];
Hlist = Hdata.Hlist%;[Hdata.reflist, Hdata.testlist];
nameslist = Hdata.dfnames;%[Hdata.refnames, Hdata.testnames];
currgmm = gmmlist{sampleInd};
currH = Hlist{sampleInd};
currname = nameslist{sampleInd};

modelcolor = [0.6350, 0.0780, 0.1840];
datacolor = [0, 0, 0];
mgray=[0.6,0.6,0.6];

n = size(currH, 2);
X = random(currgmm, n)';

f = length(fplotinds);
xy = nchoosek(1:f,2);
Hinds=randsample(size(currH,2),500);
Xinds=randsample(size(X,2),500);

fig1 = figure;
fig2 = figure;
for i=1:length(xy)
    
    x=xy(i,1);
    y=xy(i,2);
    subplotpos = (y-1)*f+x;
    
    figure(fig1)
    ax1=subplot(f,f,subplotpos);
    scatter(currH(fplotinds(xy(i,1)),Hinds), currH(fplotinds(xy(i,2)),Hinds),10,datacolor, '.'); 
    set(gca,'ycolor',mgray);
    set(gca,'xcolor',mgray);
    set(gca,'Xtick',[]);
    set(gca,'Ytick',[]);
    axis tight
    yl=ylim()
    xl=xlim()
    
    figure(fig2)
    ax2=subplot(f,f,subplotpos);
    scatter(X(fplotinds(xy(i,1)),Xinds), X(fplotinds(xy(i,2)),Xinds),10, modelcolor,'.');
    set(gca,'ycolor',mgray);
    set(gca,'xcolor',mgray);
    set(gca,'Xtick',[]);
    set(gca,'Ytick',[]);
    set(gca, 'YLim', yl);
    set(gca, 'Xlim', xl);    
end

set(fig1,'PaperOrientation','landscape');
print(fig1, '-dpdf', ['figures/',currname, '_2dproj_data.pdf'],'-bestfit')

set(fig2,'PaperOrientation','landscape');
print(fig2, '-dpdf', ['figures/',currname, '_2dproj_model.pdf'],'-bestfit')