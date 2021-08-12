function Hdata = calcTSNE(Hdata, method)
% function Hdata = calcTSNE(Hdata)
% This function calculates tsne coordinates for Hdata and returns it in the
% Hdata variable. 
% 
% Two options for 'method' are :
% 
% 'joint': computes all t-SNE coordinates for all datasets together
%          can take a long time
%          stored in Hdata.jointtsnelist
% 'independent' (default): computes all t-SNE coordinates 
%          fast computation
%          stored in Hdata.tsnelist
% 

if nargin==1
    method='independent';
end

switch(method)
    case 'independent'
        
        numsamples = length(Hdata.Hlist);
        Ylist={};
        for i=1:numsamples
            currH = Hdata.Hlist{i};
            [Y,loss] = tsne(currH');
            Ylist{i} = Y;
        end
        Hdata.tsnelist = Ylist;
        
    case 'joint'
        
        dflist = Hdata.Hlist;
        Ylist = jointTSNE(dflist);
        Hdata.jointtsnelist = Ylist;

end