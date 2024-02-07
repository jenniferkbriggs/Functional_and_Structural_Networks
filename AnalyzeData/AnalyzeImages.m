%function [T,Tpp,Tpm,Tmm,S2, S3, S4, S5, Lrand2, Erand2, Crand2, S_L2] = AnalyzeImages(k)
%Ths code was used to analyze all run network analysis on imaging data
%Jennifer Briggs, 2021

clear all
close all
%clc
load('/Volumes/Briggs_10TB/SizeDependence/Analysis/Images/VidInfo.mat');

%folderpaths = {'1_sec\', '6_sec\','10_sec\'}
%folderpaths = {'6_sec\'}
%imagepaths(1).name = {'39\', '41\', '35\', '43\'}
%imagepaths(2).name = {'20\', '22\', '23\', '24\'}
%imagepaths(3).name = {'11\', '4\', '1\', '7\', '10\','19\'}
addpath('~/Documents/GitHub/Simulations');
addpath('~/Documents/GitHub/UniversalCode');
savetime =  datestr(datetime('today'),'yyyymmdd');

path = ['/Volumes/Briggs_10TB/SizeDependence/Analysis/Images/TrywFPS/'];
savename = ['/Volumes/Briggs_10TB/SizeDependence/Analysis/Images/Network_SDfilter' savetime];
cti = 1;
%Images = {'39', '41','20','23', '1', '7', '10', '19';
            %'35', '43','22','24', '4', '11','14','15'}
% Images = {'1', '7', '8','10', '19','39', '41', '23', '24', '3','4', '11', '14', '15','35', '43', '20', '22'} 
Opts.figs = 1;
Opts.printasSubplot = 0;
Opts.multiseed = 0;
Opts.multiseednum = 0;
Opts.subplotdim = 320;
Opts.Subplotnum = 1;


folders = dir([path '*']);
j = 1;
cnt = 1;

if 1
for illy = 1:length(folders)

     Images = dir([path folders(illy).name '\*' '.lsm']);
%     %Set Th
    clearvars jt
    if length(Images)>0
    for kl = 1:length(Images)
        jt(kl) = contains(Images(kl).name, 'WT')
    end
    wt = find(jt == 1)
    Opts.figs = 1;

    for ji = 1:length(wt)
        load([path folders(illy).name '\' Images(wt(ji)).name '\' 'CaWaveForm.mat'])
        Thr = findoptRth(CellTC)
        [N, Adj, x, y,~,~,s] = links(CellTC, Thr, Opts,1);
        
        optThr(cnt) = Thr;
        avN(cnt)= mean(N);
        
        

        Fit(cnt)=s
      
        saveas(figure(2), ['D:\SizeDependence\Analysis\Images\Thropt\Wt' num2str(cnt) 'log.png'])
        saveas(figure(1), ['D:\SizeDependence\Analysis\Images\Thropt\Wt' num2str(cnt) '.png'])
        cnt = cnt+1;
        close all
    end
    end
end

T = table([1:cnt-1]', optThr', avN', Fit', 'VariableNames', {'Islet', 'Optimal Threshold', 'Average Connections', 'Fit to powerlaw'})
writetable(T, ['D:\SizeDependence\Analysis\Images\Thropt\Thr and Fits.xlsx' ])
Thr= mean(optThr);

else
    Thr = .90;
end

%% 


for illy = 1:length(folders)
       Images = dir([path folders(illy).name '/*' '.lsm']);

for i = 1:length(Images)
    load([path folders(illy).name '/' Images((i)).name '/' 'CaWaveForm.mat'])
    load([path folders(illy).name '/' Images((i)).name '/' 'Masksv2.mat'])
    %load([path folders(illy).name '\' Images((i)).name '\' 'CellNumber.mat'])

    for ll = 1:size(CellTC,2)
        [xp,yp] = find(CellMask == ll);
        x(ll,i) = mean(xp);
        y(ll,i) = mean(yp);
    end
    Folder(j).Network(i).data = (CellTC);
   Folder(j).Network(i).name = Images(i).name;
    Folder(j).Network(i).CellposX = x;
    Folder(j).Network(i).CellposY = y;
   % Network(i).Rij = Rij;
   % clear xp yp x y
    [col, row] = find(isnan(Folder(j).Network(i).data));
    Folder(j).Network(i).data(:,unique(row)) = [];

    Opts.Subplotnum = i;
    Opts.figs = 0;
    [Folder(j).Network(i).Rij, pval]=corr(Folder(j).Network(i).data);

    
% for ll = 1:numcells
%     [xp,yp] = find(CellMask == ll);
%     x(ll) = mean(xp);
%     y(ll) = mean(yp);
% end
    if 1
    [Folder(j).Network(i).N, Folder(j).Network(i).Adj, Folder(j).Network(i).x, ...
        Folder(j).Network(i).y,Folder(j).Network(i).pval, ...
        Folder(j).Network(i).Rij] = links(Folder(j).Network(i).data, Thr, Opts,0);
    else
        [Folder(j).Network(i).N, Folder(j).Network(i).Adj, Folder(j).Network(i).x, ...
        Folder(j).Network(i).y,Folder(j).Network(i).pval, ...
        Folder(j).Network(i).Rij] = linksXcor(Folder(j).Network(i).data, Thr, Opts,0);
        
    end
   Folder(j).Network(i).Upper = triu(Folder(j).Network(i).Rij,1);

%     adj = (Folder(j).Network(i).Rij>Thr);
%     adj = adj-diag(diag(adj));
%     disp(mean(sum(adj)))
adj = (Folder(j).Network(i).Adj);
%     x=NaN;
%     y = NaN;
    try
        kavg(j,i) = mean(sum(adj))/length(adj);%max(sum(adj));
    catch
        disp('Extending arrays')
        kavg_nonorm(j-1,i) = NaN;
        kavg(j-1,i) = NaN;
        Ravg(j-1,i) = NaN;
        Cellsize(j-1,i) = NaN;
        L(j-1,i) = NaN;
        Lrand(j-1,i) = NaN;
        Eglob(j-1,i) = NaN;
        Erand(j-1,i) = NaN;
        nopath(j-1,i) = NaN;
        dist(j-1,i) = NaN;
        Cavg(j-1,i) = NaN;
        Crand(j-1,i) = NaN'
    end
            kavg_nonorm(j,i) = mean(sum(adj));

   kavg(j,i) = mean(sum(adj))/length(adj);%/max(sum(adj));
   Ravg(j,i) = mean2(adj);
   Cellsize(j,i) = size(Folder(j).Network(i).data,2);
Opts.Clustering = 1;
Opts.PL = 1;
Opts.Random = 0;
Opts.GS = 0;
 [L(j,i),Lrand(j,i,:),Eglob(j,i),Erand(j,i,:),Eloc(j,i),nopath(j,i), dist(j,i), ...
     Cavg(j,i),Crand(j,i,:), SE(j,i,:),SL(j,i,:), ~, ~] = graphProp_multi(adj,x,y,1,1000,Opts);
 
 %clearvars x y
 Names(j,i)={Folder(j).Network(i).name};
 Foldername(j,i) = {folders(illy).name};
% figure, bar(Folder(j).Network(i).x',Folder(j).Network(i).y')
%  title([Foldername(j,i) Folder(j).Network(i).name])

%     T = array2table([Folder(j).Network(i).x',Folder(j).Network(i).y'],'Variablenames',{'Links', 'Perc'})
%     writetable(T,[path folders(illy).name '\' Images(i).name 'Histogram.txt'])      
%     close all
    Thrr(j,i) = Thr;
    if Cellsize(j,i) ~= 0
    Type(j,i) = Vidinfo(j+2).type(i);
    Mouse(j,i) = Vidinfo(j+2).Mouse(i);
    end
end

 if length(Images)>0
 cti = cti+i-1;
 j = j+1;
 end
end
% Cavg(isnan(Cavg))=0;
% Crand(isnan(Crand)) = 0; %If Crand is 0, the random graph had no clustering, but Cavg/0 = \infty so we adjust
% Crand(Crand == 0)=.01;
% S_L = (Cavg./Crand)./(L./Lrand);
% S = (Cavg./Crand)./(Erand./Eglob);
%S(isinf(S)) = 0;

%% 
bins = 20;
x = [];
% [row,col] = find(Type == 3) %Find all images of type selected %1 = mm. 2 = pm. 3 = pp
% Histx = linspace(0,100,bins);
% Hist = zeros(length(row), bins);
% figure
for i = 1:length(row)
    adj = (Folder(row(i)).Network(col(i)).Adj); %Pull adj matrix data
    K = sum(adj);                               %Calculate the degree of each node
    mk(i) = max(K)                              %Calculate max degree of network
        binwidth = [0:mk(i)/(bins-1):mk(i)]        
    for j = 1:length(adj)
        if j == 66
            dd = 1
        end
        lb = find(K(j) < binwidth);
        ub = find(K(j) > binwidth);
        if ~isempty(lb)
        lb = lb(1);
        if K(j)==0
        Hist(i,lb-1) = Hist(i,lb-1)+.75;
        Hist(i,lb) = Hist(i,lb)+.25;
        elseif lb(1) == bins
        Hist(i,lb-1) = Hist(i,lb-1)+.25;
        Hist(i,lb) = Hist(i,lb)+.75;
        else
        Hist(i,lb-1) = Hist(i,lb-1)+.25;
        Hist(i,lb) = Hist(i,lb)+.5;
        Hist(i,lb+1) = Hist(i,lb+1)+.25;
        end
        else
        ub = ub(end)
        Hist(i,ub) = Hist(i,ub)+.25;
        Hist(i,ub+1) = Hist(i,ub+1)+.75;
        end
        if sum(Hist(i,:)) ~= j
         dd = 1
        end
    end
%     Hist(i,:) = Hist(i,:)./length(adj);
%     hold on
%     bar(Hist(i,:))
%     drawnow
    %%
%    Folder(row(i)).Network(col(i)).name

   % x=[x sum(adj)];
end
%% 
x = (x./max(x)).*100;
Lrand2=struct
Crand2=struct
Erand2 = struct
SLst = struct
SEst = struct



for lt = [1:3]
    for jt = 1:max(unique(Mouse(Type == lt)))
    mm = find(Type == lt & Mouse == jt);
    [rm, cm] = find(Type == lt & Mouse == jt);

    loc2{lt,jt} = num2str(reshape((mm),1,[]));
    kavg2(lt,jt) = nonzeros(mean(nonzeros(kavg(mm))));
    Ravg2(lt,jt) = mean(nonzeros(Ravg(mm)));
    Cellsize2(lt,jt) = mean(nonzeros(Cellsize(mm)));
    L2(lt,jt) = mean(nonzeros(L(mm)));
    Eglob2(lt,jt) = mean(nonzeros(Eglob(mm)));
    nopath2(lt,jt) = mean(nonzeros(nopath(mm)));
    Cavg2(lt,jt) = mean((Cavg(mm)));
%    S2(lt,jt) =mean((S(mm)));
    dist2(lt,jt) = mean(nonzeros(dist(mm)));
    Eloc2(lt,jt) = mean(Eloc(mm));

        
 if ~isempty(rm)
    Lrand2(lt).Mouse(jt).data = reshape(squeeze(Lrand(rm(1),cm,:)),...
        [1,size(squeeze(Lrand(rm(1),cm,:)),2)*...
        size(squeeze(Lrand(rm(1),cm,:)),1)]);
    Erand2(lt).Mouse(jt).data = reshape(squeeze(Erand(rm(1),cm,:)),...
        [1,size(squeeze(Erand(rm(1),cm,:)),2)*...
        size(squeeze(Erand(rm(1),cm,:)),1)]);
    Crand2(lt).Mouse(jt).data = reshape(squeeze(Crand(rm(1),cm,:)),...
        [1,size(squeeze(Crand(rm(1),cm,:)),2)*...
        size(squeeze(Crand(rm(1),cm,:)),1)]);
    
    SL2(lt,jt) = mean(mean(squeeze(SL(rm(1),cm,:)),2));
    SE2(lt,jt) = mean(mean(squeeze(SE(rm(1),cm,:)),2));

    SLst(lt).Mouse(jt).data = (reshape(squeeze(SL(rm(1),cm,:)),...
        [1,size(squeeze(SL(rm(1),cm,:)),2)*...
        size(squeeze(SL(rm(1),cm,:)),1)]));
    SEst(lt).Mouse(jt).data = (reshape(squeeze(SE(rm(1),cm,:)),...
        [1,size(squeeze(SE(rm(1),cm,:)),2)*...
        size(squeeze(SE(rm(1),cm,:)),1)]));
 else
 end
end
end


save('/Volumes/Briggs_2TB/SizeDependence/Analysis/Images/Network_Aug32021.mat')
%% 
kavg2(kavg2 == 0) = NaN;
Cellsize2(Cellsize2 == 0) = NaN;
Eglob2(Eglob2 == 0) = NaN;
%L2(L2 == 0) = NaN;
nopath2(nopath2 == 0) = NaN;
Cavg2(Cavg2 == 0) = NaN;
dist2(dist2 == 0) = NaN;
%S2(S2 == 0) = NaN;
Ravg2(Ravg2 == 0) =NaN;

% T = array2table([ kavg2, Ravg2, Cellsize2, L2, Eglob2, ...
%    nopath2, Cavg2, S2, dist2]);
% writetable(T, [savename '_ByMouse.txt'])
%% MM
for lt = 1:max((unique(Mouse(Type == 1))))
        mm = find(Type == 1 & Mouse == lt);
        [rm, cm] = find(Type == 1 & Mouse == lt);

for jt = 1:length(mm)
    kavg3(lt,jt) = kavg(mm(jt));
    Foldername3(lt,jt) = cellstr(strjoin([Foldername(mm(jt)) Names(mm(jt))]));
    Ravg3(lt,jt) = (nonzeros(Ravg(mm(jt))));
    Cellsize3(lt,jt) = (nonzeros(Cellsize(mm(jt))));
    L3(lt,jt) = (nonzeros(L(mm(jt))));
    Eglob3(lt,jt) = (nonzeros(Eglob(mm(jt))));
   % nopath3(lt,jt) = (nonzeros(nopath(mm(jt))));
    Cavg3(lt,jt) = ((Cavg(mm(jt))));
    dist3(lt,jt) = (nonzeros(dist(mm(jt))));
    Eloc3(lt,jt) = mean(Eloc(mm));
    try
    S3(lt,jt) =((S(mm(jt))));
    catch
    end
end
    Lrand3(lt).data = reshape(squeeze(Lrand(rm(1),cm,:)),1,[])';
    Erand3(lt).data = reshape(squeeze(Erand(rm(1),cm,:)),1,[])';
    Crand3(lt).data = reshape(squeeze(Crand(rm(1),cm,:)),1,[])';
    
    SL3(lt,jt) = mean(mean(squeeze(SL(rm(1),cm,:)),2));
    SE3(lt,jt) = mean(mean(squeeze(SE(rm(1),cm,:)),2));

    SLst3(lt).Mouse(jt).data = (reshape(squeeze(SL(rm(1),cm,:)),...
        [1,size(squeeze(SL(rm(1),cm,:)),2)*...
        size(squeeze(SL(rm(1),cm,:)),1)]));
    SEst3(lt).Mouse(jt).data = (reshape(squeeze(SE(rm(1),cm,:)),...
        [1,size(squeeze(SE(rm(1),cm,:)),2)*...
        size(squeeze(SE(rm(1),cm,:)),1)]));


end

kavg3(kavg3 == 0) = NaN;
Cellsize3(Cellsize3 == 0) = NaN;
Eglob3(Eglob3 == 0) = NaN;
L3(L3 == 0) = NaN;
%nopath3(nopath3 == 0) = NaN;
%Cavg3(Cavg3 == 0) = NaN;
dist3(dist3 == 0) = NaN;
%S3(S3 == 0)= NaN;
Ravg3(Ravg3==0) = NaN;
% Tmm = array2table([ kavg3, Ravg3, Cellsize3, L3, Eglob3, ...
%     nopath3, Cavg3, S3, dist3]);
% writetable(Tmm, [savename '_mm.txt'])
% 
%% PM
for lt = 1:max((unique(Mouse(Type == 2))))
        mm = find(Type == 2 & Mouse == lt);
         [rm, cm] = find(Type == 2 & Mouse == lt);

for jt = 1:length(mm)
    kavg4(lt,jt) = kavg(mm(jt));
    Foldername4(lt,jt) = cellstr(strjoin([Foldername(mm(jt)) Names(mm(jt))]));
    Ravg4(lt,jt) = (nonzeros(Ravg(mm(jt))));
    Cellsize4(lt,jt) = (nonzeros(Cellsize(mm(jt))));
    L4(lt,jt) = (nonzeros(L(mm(jt))));
    Eglob4(lt,jt) = (nonzeros(Eglob(mm(jt))));

    nopath4(lt,jt) = (nonzeros(nopath(mm(jt))));
    Cavg4(lt,jt) = ((Cavg(mm(jt))))
    dist4(lt,jt) = (nonzeros(dist(mm(jt))));
    Eloc4(lt,jt) = mean(Eloc(mm));
    try
    S4(lt,jt) =((S(mm(jt))));
    catch
    end
end
if ~isempty(rm)

    Lrand4(lt).data = reshape(squeeze(Lrand(rm(1),cm,:)),1,[])';
    Erand4(lt).data = reshape(squeeze(Erand(rm(1),cm,:)),1,[])';
    Crand4(lt).data = reshape(squeeze(Crand(rm(1),cm,:)),1,[])';
    
    SL3(lt,jt) = mean(mean(squeeze(SL(rm(1),cm,:)),2));
    SE3(lt,jt) = mean(mean(squeeze(SE(rm(1),cm,:)),2));
    SLst3(lt).Mouse(jt).data = (reshape(squeeze(SL(rm(1),cm,:)),...
        [1,size(squeeze(SL(rm(1),cm,:)),2)*...
        size(squeeze(SL(rm(1),cm,:)),1)]));
    SEst3(lt).Mouse(jt).data = (reshape(squeeze(SE(rm(1),cm,:)),...
        [1,size(squeeze(SE(rm(1),cm,:)),2)*...
        size(squeeze(SE(rm(1),cm,:)),1)]));
end
end
kavg4(kavg4 == 0) = NaN;
Cellsize4(Cellsize4 == 0) = NaN;
Eglob4(Eglob4 == 0) = NaN;
L4(L4 == 0) = NaN;
nopath4(nopath4 == 0) = NaN;
%Cavg4(Cavg4 == 0) = NaN;
dist4(dist4 == 0) = NaN;
%S4(S4==0) = NaN;
Ravg4(Ravg4==0) = NaN;

% aa = [ kavg4, Ravg4, Cellsize4, L4, Eglob4, ...
%     nopath4, Cavg4, S4, dist4];
% aa(5,:) = [];
% Tpm = array2table(aa);
% writetable(Tpm, [savename '_pm.txt'])

%% PP
for lt = 1:max((unique(Mouse(Type == 3))))
        mm = find(Type == 3 & Mouse == lt);
         [rm, cm] = find(Type == 3 & Mouse == lt);

for jt = 1:length(mm)
    kavg5(lt,jt) = kavg(mm(jt));
    Foldername5(lt,jt) = cellstr(strjoin([Foldername(mm(jt)) Names(mm(jt))]));
    Ravg5(lt,jt) = (nonzeros(Ravg(mm(jt))));
    Cellsize5(lt,jt) = (nonzeros(Cellsize(mm(jt))));
    L5(lt,jt) = (nonzeros(L(mm(jt))));
    Eglob5(lt,jt) = (nonzeros(Eglob(mm(jt))));
   % nopath5(lt,jt) = (nonzeros(nopath(mm(jt))));
    Cavg5(lt,jt) = ((Cavg(mm(jt))));

    dist5(lt,jt) = (nonzeros(dist(mm(jt))));
    Eloc5(lt,jt) = mean(Eloc(mm));
    try
    S5(lt,jt) =(nonzeros(S(mm(jt))));
    catch
    end
end
if ~isempty(rm)

    Lrand5(lt).data = reshape(squeeze(Lrand(rm(1),cm,:)),1,[])';
    Erand5(lt).data = reshape(squeeze(Erand(rm(1),cm,:)),1,[])';
    Crand5(lt).data = reshape(squeeze(Crand(rm(1),cm,:)),1,[])';
    
    SL5(lt,jt) = mean(mean(squeeze(SL(rm(1),cm,:)),2));
    SE5(lt,jt) = mean(mean(squeeze(SE(rm(1),cm,:)),2));
    SLst5(lt).Mouse(jt).data = (reshape(squeeze(SL(rm(1),cm,:)),...
        [1,size(squeeze(SL(rm(1),cm,:)),2)*...
        size(squeeze(SL(rm(1),cm,:)),1)]));
    SEst5(lt).Mouse(jt).data = (reshape(squeeze(SE(rm(1),cm,:)),...
        [1,size(squeeze(SE(rm(1),cm,:)),2)*...
        size(squeeze(SE(rm(1),cm,:)),1)]));
end
end
kavg5(kavg5 == 0) = NaN;
Cellsize5(Cellsize5 == 0) = NaN;
Eglob5(Eglob5 == 0) = NaN;
L5(L5 == 0) = NaN;
nopath5(nopath5 == 0) = NaN;
%Cavg5(Cavg5 == 0) = NaN;
dist5(dist5 == 0) = NaN;
Ravg5(Ravg5==0) = NaN;

% Tpp = array2table([ kavg5, Ravg5, Cellsize5, L5, Eglob5, ...
%     nopath5, Cavg5, S5, dist5]);
% writetable(Tpp, [savename '_pp.txt'])
% % 



% Foldername2= reshape(Foldername', size(Foldername,1)*size(Foldername,2),1)
% Names2 = reshape(Names', size(Foldername,1)*size(Foldername,2), 1);
% kavg2 = num2cell(reshape(kavg', size(Foldername,1)*size(Foldername,2), 1));
% dist = num2cell(reshape(dist', size(Foldername,1)*size(Foldername,2), 1));
% 
% Ravg2 = num2cell(reshape(Ravg', size(Foldername,1)*size(Foldername,2), 1));
% Cellsize2 = num2cell(reshape(Cellsize', size(Foldername,1)*size(Foldername,2), 1));
% L2 = num2cell(reshape(L', size(Foldername,1)*size(Foldername,2), 1));
% Lrand2= num2cell(reshape(Lrand', size(Foldername,1)*size(Foldername,2), 1));
% Eglob2 = num2cell(reshape(Eglob', size(Foldername,1)*size(Foldername,2), 1));
% Erand2 = num2cell(reshape(Erand', size(Foldername,1)*size(Foldername,2), 1));
% nopath2 = num2cell(reshape(nopath', size(Foldername,1)*size(Foldername,2), 1));
% Cavg2 = num2cell(reshape(Cavg', size(Foldername,1)*size(Foldername,2), 1));
% Crand2 = num2cell(reshape(Crand', size(Foldername,1)*size(Foldername,2), 1));
%  Thrr2 = num2cell(reshape(Thrr', size(Foldername,1)*size(Foldername,2), 1));
%   kavg_nonorm2 = num2cell(reshape(kavg_nonorm', size(Foldername,1)*size(Foldername,2), 1));
% S2 = num2cell(reshape(S', size(Foldername,1)*size(Foldername,2), 1));
% T = array2table([Foldername2,Names2, kavg2, Ravg2, Cellsize2, L2, Eglob2, ...
%     nopath2, Cavg2, S2, Thrr2, kavg_nonorm2,dist],...
%     'VariableNames', ...
%     {'Folder', 'Islet', 'Connections', 'Ravg', 'Cellsize' ,'Path Length', 'Efficiency',...
%     'Cells w/no path', 'Clustering', 'SmallWorldness', 'Threshold','Kavg not normalized', 'Average distance btw connected cells'})
%  writetable(T,[savename '.xlsx'])      
%end
save('MuiltipAug9.mat')
