%Network Analysis for Size Dependence%
clear all
close all
mac = 1
if mac == 1
addpath('~/Documents/GitHub/UniversalCode/');
addpath('~/Documents/GitHub/Simulation-Analysis/');
addpath('~/Documents/GitHub/Simulations/');
addpath('~/Documents/GitHub/Network-Analysis/');

%%%%%%%%%%%%%%%%%%%%%%%%%%
savetime =  datestr(datetime('today'),'yyyymmdd');
freqcalc = 0
%datapath = 'D:\Islet\SizeDependance\Size_depend'
%dataset2 = {'50cell','100cell','300cell','1000cell','2000cell','3000cell'}
%dataset2 = {'FullCoup/1000cell', 'HalfCoup_0311/1000cell', 'QuartCoup_0311/1000cell'}
%dataset2 = {'FullCoup\1000cell', 'FullCoup\1000cell', 'FullCoup\1000cell', 'FullCoup\1000cell', 'FullCoup\1000cell'}%,'HalfCoup_0311\1000cell','QuartCoup_0311\1000cell'}
dataset2 = {'FullCoup/1000cell'}
dataname = {'/seed1/', '/seed2/','/seed3/','/seed4/','/seed5/'}
datapath = '/Volumes/Briggs_2TB/SizeDependence'
slash = '/';
else
addpath('~/UniversalCode/');
slash = '/';
%addpath('~/Simulation-Analysis/');
addpath('~/Simulations/');
%%%%%%%%%%%%%%%%%%%%%%%%%%
savetime =  datestr(datetime('today'),'yyyymmdd');
freqcalc = 0
%datapath = 'D:\Islet\SizeDependance\Size_depend'
%dataset2 = {'50cell','100cell','300cell','1000cell','2000cell','3000cell'}
dataset2 = {'FullCoup/1000cell', 'FullCoup/1000cell', 'FullCoup/1000cell'}%, 'FullCoup/1000cell', 'FullCoup/1000cell'}%,'HalfCoup_0311\1000cell','QuartCoup_0311\1000cell'}
dataname = {'/seed1/', '/seed2/','/seed3/','/seed4/','/seed5/'}
datapath = '~/Documents/SizeDependenceIslet/'  
end
%Names = {'All', 'P2', 'P1'}
%Names = {'.9997', '.9998', '.9999', '.99991' '.99992'}
Names = {'gCoup = .12 pS' ,'gCoup = .06 pS', 'gCoup = .03 pS'};
%Names = {'P2'}
Threshold = .9995; %Set R_th for network analysis
%Threshold = [.9997, .9998, .9999, .99991, .99992];
savename = ['2nd_PhaseXcor' savetime];
%Set options for network analysis
Opts.figs = 0;
Opts.change_thres = 0;
Opts.chage_phase = 0;
Opts.printasSubplot = 1;
%Set the last number to 0 always, first number is rows second is columns
Opts.subplotdim = 310;
Opts.multiseed = 1;

fig = 1
seeds = length(dataname);
Isizes = length(dataset2);
%% Import Data
% try 
%     for l = 1:length(dataset2)
%     for i = 1:length(dataname)
%     load([savename '.mat'])
%     end
%     end
% catch
ca = struct;
Rvars = struct;
tic
for l = 1:length(dataset2)
    for i = 1:length(dataname)
    disp(['Importing ' cell2mat(dataset2(l))])
    capath = fullfile(datapath, dataset2(l), dataname(i), 'calcium.txt');
    casize = importdata(cell2mat(capath));
    if i == 1
       ca2(l).data(:,:,i) = casize;
    else
       sizes = [size(casize,1), size(ca2(l).data(:,:,i-1),1)]
       [sizemin,seedmin] = min(sizes)
       for tt = 1:l
       ca2(tt).data(sizemin+1:end,:,:)=[]
       end
       casize(sizemin+1:end,:)=[];
       ca2(l).data(:,:,i) = casize;
    end
    ca2(l).name = dataset2(l);
    
    pospath = fullfile(datapath, dataset2(l), dataname(i), 'XYZpos.txt');
    pos(l).data(:,:,i) = importdata(cell2mat(pospath));
    pos(l).name = dataset2(l);
    
    connpath = fullfile(datapath, dataset2(l), dataname(i), 'NN10A.txt');
    conn(l).data(:,:,i) = importdata(cell2mat(connpath));    
    calcium = ca2(l).data(:,:,i);    
    ca(l).data(:,:,i) = calcium(2000:end,:);
%     if l == 1
% %     calcium = importdata(cell2mat(connpath));    
%     ca(l).data(:,:,i) = calcium;
%     elseif l == 2
%     ca(l).data(:,:,i) = calcium(2000:end,:);
%     elseif l == 3
%     ca(l).data(:,:,i) = calcium(1:750,:);
%     end
    varspath = fullfile(datapath, dataset2(l), dataname(i), 'RandomVars.txt');
    Rvars(l).data(:,:,i) = importdata(cell2mat(varspath));
    Rvars(l).name = dataset2(l);
    fig1 = figure(1)%, subplot(Opts.subplotdim+l)
    
    if i == 1
    subplot(Opts.subplotdim + l)
    plot(ca(l).data), xlabel('Time'), ylabel('Ca'), title(dataset2(l))
    end
    
    

    
    end
end
toc
% try
%     ca(2).data(4997:end,:,:)=[];
% catch
%     disp('error')
% end
clear casize sizes sizemin seed min tt

fig1.Position = [325 32 1132 892.5000];
% 
%% Run Network Analysis
Links2 = struct;
Links2.N = struct;
Links2.sort = struct;
Links2.cellsort = struct;
fig = fig+1;
for l = 1:length(dataset2)
    for i = 1:length(dataname)
     Opts.printasSubplot = 1;
        Opts.Subplotnum = l;
    
    Opts.multiseednum = i;
    [Links2.N(l).data(:,i), adjm(l).data(:,:,i), kperc(l).Islet(i).data , histArrayPercShort(i).data, Rij(l).Islet(i).data(:,:) ] = links(ca(l).data(:,:,i), Threshold,Opts,fig); %%This is where network analysis is perfromed
%        [Links2.N(l).data(:,i), adjm(l).data(:,:,i), kperc(l).Islet(i).data , histArrayPercShort(i).data] = linksXcor(ca(l).data(:,:,i), Threshold,Opts,fig,cormat); %%This is where network analysis is perfromed

   [Links2.sort(l).data(:,i), Links2.cellsort(l).data(:,i)]= sort(Links2.N(l).data(:,i));
    end
end

% end

%% Calculate phase
phaseinfo = struct
phaseinfo.cell = struct;
if 0
for i = 1:Isizes
    for l = 1:seeds
%        disp(['Calculating ' cell2mat(dataset2(i))])
       [ phaseinfo.phase(i).data(:,l),cormat(:,:,l),cormattime(:,:,l)] = phase(ca(i).data(:,:,l),0,0); %time is in ms
        [psort,pcell]  = sort(phaseinfo.phase(i).data(:,l));
        phaseinfo.sort(i).data(:,l) = psort;
        phaseinfo.cell(i).data(:,l) = pcell;

    end
end
end

%% Plot for Kglyc
cellperc = .1;
Kglyc = struct;
Kglyc.sort = struct;
Kglyc.cell = struct;
Kglyc.kglyc = struct;
Kglyc.high = struct;
Kglyc.rand = struct;
Kglyc.low = struct;

fig = fig+2
for l = 1:length(dataset2)
    for i = 1:length(dataname)
    [Kglyc.sort(l).data(:,i), Kglyc.cell(l).data(:,i)]=sort(Rvars(l).data(:,10,i)); %Sorts cells in ascending order based on Kglyc
    Kglyc.kglyc(l).data(:,i) = (Rvars(l).data(:,10,i));
    [cellnum] = length(Rvars(l).data(:,10,i));
    numcell = cellperc*cellnum;
    randomcell = randi(cellnum,numcell,1);
    N = Links2.N(l).data(:,i);
    cell = Kglyc.cell(l).data(:,i);
    %%%%%%%%%%%
    Kglyc.high(l).data = mean(N(cell(end-numcell:end)));
    Kglyc.rand(l).data = mean(N(cell(randomcell)));
    Kglyc.low(l).data = mean(N(cell(1:numcell)));

    Links2.high(l).data(:,i) = (Links2.cellsort(l).data(end-numcell:end,i));
    Links2.rand(l).data(:,i) = (Links2.cellsort(l).data(randomcell,i));
    Links2.low(l).data(:,i) = (Links2.cellsort(l).data(1:numcell,i));
    a = find(kperc(l).Islet(i).data  > 60);
    hubthreshold = (a(1));
    b = find(kperc(l).Islet(i).data  <= 60);
    connect = sum(adjm(l).data(:,:,i));
    Hubs(l).hubby(i).data = find(connect>hubthreshold);%Number of cells linked to more than 50% of Islet
    hubby = Hubs(l).hubby(i).data;
    %     G = graph(adjm(l).data(:,:,i), 'omitselfloops');
%     c = centrality(G, 'betweenness');
%     cnorm = c/max(c);
%     fig = fig+1
%     figure(fig)
%     histogram(cnorm)
%     hubby = find(cnorm > .6);
%     hublow = find(cnorm <= .6);
%    	Hubs(l).hubby(i).data=hubby;
% 
%     sam = intersect(hubby, find(connect>hubthreshold))
%     title(['Seed ' num2str(i) ' Distribution of centrality'])
%     text(.5, 400,['Percent of Sustained Hubs = ' num2str(length(sam)/length(hubby)) ' & ' num2str(length(sam)/length(Hubs(l).hubby(i).data)) 'from original'])
	hublow = find(connect<=hubthreshold);
    Hubs(l).hublow(i).data = hublow%find(connect<b(end));%Number of cells linked to more than 50% of Islet
    hubmid = [1:1000];%Number of cells linked to more than 50% of Islet
    hubmid(hublow) = [];
    for iii = 1:length(hubby)
    hubmid(find(hubmid == hubby(iii))) = [];
    end
    Hubs(l).hubmid(i).data = hubmid;

    Links2.high(l).data(:,i) = (Links2.cellsort(l).data(end-numcell:end,i));
    Links2.rand(l).data(:,i) = (Links2.cellsort(l).data(randomcell,i));
    Links2.low(l).data(:,i) = (Links2.cellsort(l).data(1:numcell,i));
    
    Kglyc.hubs(l).data(:,i) = mean(Kglyc.kglyc(l).data(hubby,i));
%     Kglyc.hubsmid(l).data(:,i) = mean(Kglyc.kglyc(l).data(hubmid,i));
    Kglyc.hubslow(l).data(:,i) = mean(Kglyc.kglyc(l).data(hublow,i));
        Conny = sum((conn(1).data(:,:,i)~=-1)');
   TotConn.hubs(l).data(:,i) =mean(Conny(hubby));
%     Kglyc.hubsmid(l).data(:,i) = mean(Kglyc.kglyc(l).data(hubmid,i));
    TotConn.hubslow(l).data(:,i) = mean(Conny(hublow));
    

    Katp.hubs(l).data(:,i) = mean(Conny(hubby));
    Katp.hubslow(l).data(:,i) = mean(Conny(hublow));

    Kglyc.high(l).data(:,i) = mean(Kglyc.kglyc(l).data(Links2.high(l).data,i));
    Kglyc.rand(l).data(:,i) = mean(Kglyc.kglyc(l).data(Links2.rand(l).data,i));
    Kglyc.low(l).data(:,i) = mean(Kglyc.kglyc(l).data(Links2.low(l).data,i));
   
    gjconn = conn(1).data(:,:,i);
    gjconn = gjconn + 1;
    
    for t = 1:size(gjconn,1)
        connections = gjconn(t,:);
        connections(connections == 0) = [];
        for tt = 1:length(connections)
            conduct(tt) = Rvars(l).data(connections(tt),2,i);
        end          
        GJconduct(i).data(t) = sum(0.5*(conduct + Rvars(l).data(t,2,i)));
                clear conduct

    end
     
    Gcoup.hubs(l).data(:,i) = mean(GJconduct(i).data(hubby));
%     Gcoup.hubsmid(l).data(:,i) = mean(Rvars(l).data(hubmid,2,i));
    Gcoup.hubslow(l).data(:,i) = mean(GJconduct(i).data(hublow));

    Gcoup.high(l).data(:,i) = mean(Rvars(l).data(Links2.high(l).data,2,i));
    Gcoup.rand(l).data(:,i) = mean(Rvars(l).data(Links2.rand(l).data,2,i));
    Gcoup.low(l).data(:,i) = mean(Rvars(l).data(Links2.low(l).data,2,i));
    end
    %%%%%%%%%%%  
    if Opts.figs
fig = fig+1
figure(fig)
%     if length(dataset2) >1
%     subplot(Opts.subplotdim+l)
%     end
    [h] = multiseedplotsave(fig,[Gcoup.high(l).data; Gcoup.rand(l).data; Gcoup.low(l).data], ...
        ['Average Coupling Conductance for cells sorted by links_' Names{l}], ...
        {'Highest 10% Links', 'Rand 10% Links', 'Lowest 10% Links'},...
        {'gCoup [pS]'}, [savename '_Gcoup' Names{l}])
%     set(gca,'XTickLabel',{'Hubs (>60% synchronization)', 'Mid Links (30% - 60%)', 'Low Links (< 30% synchronization)'});
%     title({'Average Coupling Conductance for cells sorted by links' '(Islet Size = ' num2str(cellnum) 'cells)'})
%     ylabel({'gCoup [pS]'})
 fig = fig +1;
[h] = multiseedplotsave(fig,[Gcoup.hubs(l).data;Gcoup.hubslow(l).data], ...
        ['Katp' Names{l}], ...
        {'Hubs (>60% synchronization)', 'Low Links (< 20% synchronization)'},...
        {'Katp '}, [savename 'Katp' Names{l}])
    fig = fig +1;
    
    fig = fig +1;
[h] = multiseedplotsave(fig,[Katp.hubs(l).data;Katp.hubslow(l).data], ...
        ['Average Coupling Conductance for cells sorted by links' Names{l}], ...
        {'Hubs (>60% synchronization)', 'Low Links (< 20% synchronization)'},...
        {'gCoup [pS]'}, [savename '_Gcouphubs' Names{l}])
    fig = fig +1;
figure(fig)
%     if length(dataset2) >1
%     subplot(Opts.subplotdim+l)
%     end
    [h] = multiseedplotsave(fig,[Kglyc.hubs(l).data; Kglyc.hubslow(l).data],...
        ['Average Kglyc for cells sorted by links' Names{l}], {'Hubs (>60% synchronization)',  'Low Links (< 20% synchronization)'},...
        'Kglyc', [savename '_Kglychubs' Names{l}]);
        fig = fig +1;
figure(fig)
    [h] = multiseedplotsave(fig,[Kglyc.high(l).data; Kglyc.rand(l).data; Kglyc.low(l).data],...
        ['Average Kglyc for cells sorted by links' Names{l}], {'Highest 10% Links', 'Rand 10% Links', 'Lowest 10% Links'},...
        'Kglyc', [savename '_Kglyc' Names{l}]);
%     set(gca,'XTickLabel',{'Hubs (>60% synchronization)', 'Mid Links (20% - 60%)', 'Low Links (< 20% synchronization)'});
%     title({'Average Kglyc for cells sorted by links' '(Islet Size = ' num2str(cellnum) 'cells)'})
%     ylabel({'Kglyc'})
    end
    clear N cell
%     xlab = {'Highest 10% Links', 'Rand 10% Links', 'Lowest 10% Links'}
end

% 
% hubby = Hubs(l).hubby(i).data;
% lowwi = Hubs(l).hublow(i).data;
% % lowwi = Links2.low(l).data(:,i);
% calcium = ca(l).data(:,:,i);
% time = [1:length(calcium)]*100/1000;
% callow = [mean(calcium(:,lowwi(1:length(lowwi)/4)),2),...
%     mean(calcium(:,lowwi(length(lowwi)/4:2*length(lowwi)/4)),2),...
%     mean(calcium(:,lowwi(2*length(lowwi)/4:3*length(lowwi)/4)),2),...
%     mean(calcium(:,lowwi(3*length(lowwi)/4:end)),2)];
% % calrand = [mean(calcium(:,randii(1:length(randii)/2)),2), mean(calcium(:,randii(length(randii)/2+1:end)),2)]
% calhub = [mean(calcium(:,hubby(1:length(hubby)/2)),2), mean(calcium(:,hubby(length(hubby)/2+1:end)),2)]
% data = array2table([time', calhub,callow], 'VariableNames', {'Time','Hub1','Hub2','Low1','Low22','Low3','Low4'})
% writetable(data, [pwd '\' savename '_' savetime '.xlsx'])
% figure
% plot(time, mean(calhub,2), 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth',2)
% hold on,plot(time, callow,'Color', [.175 .54 .60])
% ylabel('Calcium [\muM]')
% title('Calcium Dynamics')
% xlabel('Time (seconds)')
% legend('Hub Cell', 'Non-hub cell')
% 
% saveas(gcf, [pwd '\' savename '_' savetime '.fig'])
% saveas(gcf, [pwd '\' savename '_' savetime '.png'])

% %% Duty Cycle for hubs
% for l = 1:Isizes
%     for i = 1:seeds
%         [cellnum] = length(Rvars(l).data(:,10));
%         threshholdDC = 3.4759e-04
%         [dutycycle_WT(l).data(:,i), amt_on(l).data(:,i)] = getDC(ca(l).data(:,:,i), threshholdDC);
%         %amt_on(l).data = normalize(amt_on(l).data)
%         Active.high(l).data(:,i) = mean(amt_on(l).data(Links2.high(l).data,i));
%         Active.rand(l).data(:,i) = mean(amt_on(l).data(Links2.rand(l).data,i));
%         Active.low(l).data(:,i) = mean(amt_on(l).data(Links2.low(l).data,i));
%         
%         totamtactive(l,i) = mean(amt_on(l).data(:,i));
%     end
% end
% 
% 
% fig = fig+1
% 
% if length(dataset2) == 2
% figure(fig)
% h = multiseedplot2(fig,totamtactive);
% set(gca,'XTickLabel',Names);
% title({'Amount of time active' 'Threshold = ' threshholdDC})
% ylabel({'Amount Active'})
% elseif length(dataset2)
% figure(fig)
% h = multiseedplot6(fig,totamtactive);
% set(gca,'XTickLabel',Names);
% title({'Amount of time active' 'Threshold = ' threshholdDC})
% ylabel({'Amount Active'})
% end
%% Probability that a cell is linked given functional connection
clear t counter
for i = 1:Isizes
   for l = 1:seeds
       tic
%  for i = 4
%      for l = 2
% hubby = Hubs(i).hubby(l).data;
posit = pos(i).data(:,:,l);
if length(Threshold) == 1
        gjconn = conn(i).data(:,:,l);
        gjconn = gjconn + 1; % Data is in CPP so we need to start with 1 instead of 0
else
    gjconn = conn(1).data(:,:,l);
    gjconn = gjconn + 1;
end
        avgjconn(i,l) = sum(nnz(gjconn))/size(gjconn,1);
        adj = adjm(i).data(:,:,l);
        adj = adj - diag(diag(adj)); % Make adj matrix diags 0

  %Convert gjconn to look like adj
  gj2adj = zeros(size(adj));
t = 1;
if 1 %Create a network of gap junction connected cells
    for c = 1:size(gjconn,1)
        for n=1:size(gjconn,2)
%             disp([c n])
%             disp(gjconn(c,n))
            if gjconn(c,n) ~= 0
                 gj2adj(c,gjconn(c,n)) = 1;
                 counter(t) = adj(c,gjconn(c,n));
                 if isnan(counter(t))
                     disp('oops')
                 end
                 t=t+1;
            end
        end  
    end
end
posit = pos(1).data(:,:,l);
if 1
dist = sqrt((posit(:,1)-posit(:,1)').^2+(posit(:,2)-posit(:,2)').^2+...
(posit(:,3)-posit(:,3)').^2);
dist(dist==0)=100;
end
target = mean2(Rvars(i).data(:,10,l)) +1*std(Rvars(i).data(:,10,l))
if 0 %Create a network of kglyc cells > some value
    disty = []

    for c = 1:1000
            connectq = gjconn(c,:);
            connectq = unique(nonzeros(connectq));
            connectq1 = gjconn(connectq,:);
            connectq1 = unique(nonzeros(connectq1));
            connectq2 = gjconn(connectq1,:);
            connectq2 = unique(nonzeros(connectq2));
            connectqall = unique([connectq; connectq1; connectq2]);
            connectqall(find(connectqall == c))=[];
            counter = [];
            disty = [disty; dist(c,connectqall)'];
         for n=1:1000
%             if nonzeros(connectqall == n)
%                counter = [counter n];
            if 0
            if Rvars(i).data(c,10,l) > target| Rvars(i).data(n,10,l) > target 
            gj2adj(c,n) = 1;%.5*(Rvars(i).data(c,2,l)+Rvars(i).data(n,2,l));         
            end
            else %if abs((Rvars(i).data(c,10,l)-Rvars(i).data(n,10,l)))<3.6056e-05 & n~=c&dist(c,n)<.25;
            gj2adj(c,n) = abs((Rvars(i).data(c,10,l)-Rvars(i).data(n,10,l)));
            end
%              end
            
        end  
    end
 
  gj2adj(dist > .15) = 0;
  %gj2adj = (gj2adj>(mean2(nonzeros(gj2adj)) + 0.5*std(nonzeros(gj2adj),0,'all')));
  gj2adj = gj2adj<(mean2(nonzeros(gj2adj))-0*std(nonzeros(gj2adj),0,'all')) & gj2adj>0;
  gj2adj = double(gj2adj);
end
% 
% gj2adj = [1 1 0 1; 1 1 1 0; 0 1 1 0; 1 0 0 1];
% adj = [1 1 0 1; 1 1 1 0; 0 1 1 0; 1 0 0 1];
    mm(l) = mean(sum(gj2adj, 'omitnan'))
    %gjconduct = Rvars(i).data(:,2,l); %Weight by Metabolism
    gjconduct = Rvars(i).data(:,10,l);  %Weight by GJ
    [distweighted(:,:,l), adjw] = weightednet(gj2adj, gjconduct,0); %outputs graph weighted by inverse of conductance
    %[Length, syncD_all,nonsyncD_all] = WeightedDistEqualbasedonCell(adj,adjw,l)

    [col,row] = find(adj);
    dd = distweighted(:,:,l);
    dd(isinf(dd))= NaN;
    dd(dd==0) = NaN;
    for jj = 1:length(col)
    ShortestpathSync(jj) = dd(col(jj),row(jj));
    end
    
    numcells = length(adj)
%     for tgl = 1:numcells
%     GJed =  find((gj2adj(tgl,:)));
%     SYNCed = find(adj(tgl,:));
%     gjNOsync = setdiff(GJed,SYNCed);
%     syncNOgj = setdiff(SYNCed,GJed);
%     syncANDgj = intersect(SYNCed,GJed);
% 
%     for jj = 1:length(gjNOsync)
%     Shortestpath_gjNOsync(jj,tgl) = dd(tgl,gjNOsync(jj));
%     end
%     for jj = 1:length(SYNCed)
%     Shortestpath_sync(jj,tgl) = dd(tgl,SYNCed(jj));
%     end
%     for jj = 1:length(syncNOgj)
%     Shortestpath_syncNOgj(jj,tgl) = dd(tgl,syncNOgj(jj));
%     end
%      for jj = 1:length(syncANDgj)
%     Shortestpath_syncANDgj(jj,tgl) = dd(tgl,syncANDgj(jj));
%     end
%     end
%     
    x=[];
y=[];
if 0
 [Net.L(i,l),Net.Lrand(i,l), Net.Eglob(i,l), ... 
     Net.Erand(i,l), Net.Eloc(i,l), Net.nopath(i,l),...
     ~, Net.Cavg(i,l), Net.Crand(i,l)] = graphProp(adj,x,y,0)
    S_stz(i,l) = (Net.Cavg(i,l)/Net.Crand(i,l))/(Net.Erand(i,l)/Net.Eglob(i,l));
    Lsm(i,l) =  Net.L(i,l)/Net.Lrand(i,l);
    Csm(i,l) = Net.Cavg(i,l)/Net.Crand(i,l);
    Esm(i,l) = Net.Erand(i,l)/Net.Eglob(i,l);
    S_strog(i,l) = (Net.Cavg(i,l)/Net.Crand(i,l))/(Net.L(i,l)/Net.Lrand(i,l));



    shortpathsync(l) = mean(ShortestpathSync, 'omitnan');
    shortpathnosync(l) = mean(mean(dd, 'omitnan'), 'omitnan');

    mShortestpath_gjNOsync(l) = mean(mean(Shortestpath_gjNOsync, 'omitnan'), 'omitnan');
    mShortestpath_syncNOgj(l) = mean(mean(Shortestpath_syncNOgj, 'omitnan'), 'omitnan');
    mShortestpath_syncANDgj(l) = mean(mean(Shortestpath_syncANDgj, 'omitnan'), 'omitnan');

    mShortestpath_sync(l) = mean(mean(Shortestpath_sync, 'omitnan'), 'omitnan');

    clear ShortestpathSync  
    clear Shortestpath_syncANDgj Shortestpath_syncNOgj Shortestpath_gjNOsync ShortestpathSync
end
    Pr_gj(i,l) = (sum(nnz(gj2adj))/2)/(((size(gj2adj,1)-1)*(size(gj2adj,1))/2)); %Probability of gj = tot connections/tot possible connections
   %Pr_gj(i,l) = (sum(nnz(gjconn))/2)/(15*size(gjconn,1)); %Probability of gj = tot connections/tot possible connections
   Pr_syn(i,l) = (sum(nnz(adj))/2)/(((size(adj,1)-1)*(size(adj,1))/2)); %Probability of gj = tot connections/tot possible connections
%adj and gjconn are arrays where a 1 represents a connection between two
%cells. Here, we take the total number of 1's and divide by 2 (because each
%connection is represented twice) to get the total number of links in the
%Islet
    

    % calculate the graph properties given the adjacency matrix
    [Net.L(i,l),Net.EGlob(i,l),Net.CClosed(i,l),Net.ELocClosed(i,l),Net.COpen(i,l),Net.ELocOpen(i,l), Net.nopath(i,l)] = graphProperties(adj);
    
    adj2 = adj;
    gj2adj2 = gj2adj;
    
    adj(adj == 0) =NaN;
    gj2adj(gj2adj == 0) =NaN;
    
 
    both = (gj2adj == adj);
%     if sum(counter) ~= sum(nnz(both))
%         disp('Error')
%         pause(5)
%     end
    
%     hubby = Hubs(i).hubby(l).data;
%     hublow = Hubs(i).hublow(l).data;
% %     hubmid = Hubs(i).hubmid(l).data;
% 
% 
%     Pr_gjhub(i,l,1) = mean(sum(gj2adj(hubby,:),2, 'omitnan')/(((size(adj,1)-1)*(size(adj,1))/2))); %Probability of gj = tot connections/tot possible connections
% %     Pr_gjhub(i,l,2) = mean(sum(gj2adj(hubmid,:),2, 'omitnan')/(((size(adj,1)-1)*(size(adj,1))/2))); %Probability of gj = tot connections/tot possible connections
%     Pr_gjhub(i,l,2) = mean(sum(gj2adj(hublow,:),2, 'omitnan')/(((size(adj,1)-1)*(size(adj,1))/2))); %Probability of gj = tot connections/tot possible connections
% 
%     Pr_bothhub(i,l,1) = mean(sum(both(hubby,:),2, 'omitnan')/(((size(both,1)-1)*(size(both,1))/2)));
% %     Pr_bothhub(i,l,2) = mean(sum(both(hubmid,:),2, 'omitnan')/(((size(both,1)-1)*(size(both,1))/2)));
%     Pr_bothhub(i,l,2) = mean(sum(both(hublow,:),2, 'omitnan')/(((size(both,1)-1)*(size(both,1))/2)));
% 
%    Pr_synhub(i,l,1) = mean(sum(adj(hubby,:),2, 'omitnan')/(((size(adj,1)-1)*(size(adj,1)))/2)); %Probability of gj = tot connections/tot possible connections
% %    Pr_synhub(i,l,2) = mean(sum(adj(hubmid,:),2, 'omitnan')/(((size(adj,1)-1)*(size(adj,1)))/2)); %Probability of gj = tot connections/tot possible connections
%    Pr_synhub(i,l,2) = mean(sum(adj(hublow,:),2, 'omitnan')/(((size(adj,1)-1)*(size(adj,1)))/2)); %Probability of gj = tot connections/tot possible connections
% 
% 
% 
%     Pr_syn_gvn_gjhub(i,l,1) = Pr_bothhub(i,l,1)/Pr_gjhub(i,l,1);
%     Pr_gj_gvn_synhub(i,l,1) = Pr_bothhub(i,l,1)/Pr_synhub(i,l,1);
% %     Pr_syn_gvn_gjhub(i,l,3) = Pr_bothhub(i,l,3)/Pr_gjhub(i,l,3);
% %     Pr_gj_gvn_synhub(i,l,3) = Pr_bothhub(i,l,3)/Pr_synhub(i,l,3);
%     Pr_syn_gvn_gjhub(i,l,2) = Pr_bothhub(i,l,2)/Pr_gjhub(i,l,2);
%     Pr_gj_gvn_synhub(i,l,2) = Pr_bothhub(i,l,2)/Pr_synhub(i,l,2);
%         Hubcheck = 1;

    Pr_both(i,l) = (sum(nnz(both))/2)/(((size(both,1)-1)*(size(both,1))/2));
%     Pr_both2(i,l) = (sum(nnz(counter)))/(length(counter));

    Pr_syn_gvn_gj(i,l) = Pr_both(i,l)/Pr_gj(i,l)
    Pr_gj_gvn_syn(i,l) = Pr_both(i,l)/Pr_syn(i,l)
%     if length() > 1
%     [cellnum] = length(Rvars(1).data(:,10,l));
%     else
%     [cellnum] = length(Rvars(i).data(:,10,l));
%     end
    Cpos =   (((size(adj,1)-1)*(size(adj,1))/2));
    Links2.Ci(i).data(:,l) = Links2.N(i).data(:,l)./Cpos;
   % Cavg(i,l) = mean(Links2.Ci(i).data(:,l));

    clear counter
    toc
    end
end

%% Location of hubs
% icenter = [0.5, 0.5 0.5];
% for l = 1:Isizes
%     for i = 1:seeds
%      [cellnum] = length(Rvars(l).data(:,10,i));
%     distfromcoord = ((pos(l).data(:,:,i)-icenter));
%     distfromcent(l).data(:,i) = (distfromcoord(:,1).^2+distfromcoord(:,2).^2+distfromcoord(:,3).^2).^(1/2)
%     
% %     disted = ((pos(l).data-[1,1,1]));
% %     distfromedge(l).data = (disted(:,1).^2+disted(:,2).^2+disted(:,3).^2).^(1/2)
% %     
%     dist.high(l).data(:,i) = mean(distfromcent(l).data(Links2.high(l).data,i));
%     dist.rand(l).data(:,i) = mean(distfromcent(l).data(Links2.rand(l).data,i));
%     dist.low(l).data(:,i) = mean(distfromcent(l).data(Links2.low(l).data,i));
% % 
% %     diste.high(l).data = mean(distfromedge(l).data(Links2.high(l).data));
% %     diste.rand(l).data = mean(distfromedge(l).data(Links2.rand(l).data));
% %     diste.low(l).data = mean(distfromedge(l).data(Links2.low(l).data));
%     end
%     fig7 = figure(7)
%     fig7.Position = [325 32 1132 892.5000];
%     subplot(Opts.subplotdim+l)
%     [h] = multiseedplot(7,[(dist.high(l).data); dist.rand(l).data; (dist.low(l).data)])
% 
%     set(gca,'XTickLabel',{'High Links', 'Rand Links', 'Low Links'});
%     title({'Average distance from center sorted by links' '(Islet Size = ' num2str(cellnum) 'cells)'})
%     ylabel({'Distance from center'})
%     
%     fig8 = figure(8)
%     fig8.Position = [325 32 1132 892.5000];
%     subplot(Opts.subplotdim+l)
%     bar([mean(diste.high(l).data), mean(diste.rand(l).data), mean(diste.low(l).data)]);
%     set(gca,'XTickLabel',{'High Links', 'Rand Links', 'Low Links'});
%     title({'Average distance from edge sorted by links' '(Islet Size = ' num2str(cellnum) 'cells)'})
%     ylabel({'Distance from edge'})
%end


%% Extra interesting - see how this changes with thresholds
if Opts.figs
fig = fig+1
figure(fig)
KG_table_all = [mean(Links2.N(1).data)];
for iii = 2:length(dataset2)
KG_table_all = [KG_table_all; mean(Links2.N(iii).data)];
end
[h] = multiseedplot6(fig,KG_table_all)
set(gca,'XTickLabel',dataset2);
title({'Average number of links in Islet'})
ylabel({'Average number of links'})

fig = fig+1
figure(fig)
KG_table_high =[ max(Links2.N(1).data)];
for iii = 2:length(dataset2)
KG_table_high = [KG_table_high; max(Links2.N(iii).data)];
end
[h] = multiseedplot6(fig,KG_table_high)
set(gca,'XTickLabel',Names);
title('Max connections')
ylabel('Highest amount of connections in Islet')
end
%% Too dependent on size
% figure(10)
% KG_table_all = [mean(Links2.N(1).data)/50;mean(Links2.N(2).data)/100;mean(Links2.N(3).data)/300;...
%     mean(Links2.N(4).data)/1000;mean(Links2.N(5).data)/2000;mean(Links2.N(6).data)/3000];
% [h] = multiseedplot6(10,KG_table_all)
% set(gca,'XTickLabel',dataset2);
% title({'Average cell percent of connectivity'})
% ylabel({'Average number of links/total number of cells in Islet'})
% hold on, plot(mean(KG_table_all'))

%% Consider a hub a cell with 

for i = 1:Isizes
    for l = 1:seeds
        N = Links2.N(i).data(:,l);
        histArray=zeros(size(N))'; % prealocate
% a forloop to count how many times you encounter in a particular value:
            for n=1:length(N)
                histArray(1,N(n)+1)=histArray(1,N(n)+1)+1; % every time you meet the particular value, you add 1 into to corresponding bin
            end


    histArray = histArray(2:end); % removing first number corresponding to 0 links
    histArrayPerc=histArray.*100/sum(histArray); % converting into % via dividing # of cell with specific # of links by total # of links 
    m=find(histArray)    % returns all indexes of non-zero elements of the array
    maxNonZeroLinks=max(m);   % number of links to consider when normalizing probabilty
    k=1:1:maxNonZeroLinks;            % index of # of links (starting with 0 links)
    kpercent=k.*100/(maxNonZeroLinks);  % convert # of links into % of links
    histArrayPercShort = histArrayPerc(1:maxNonZeroLinks);   % cropping the hisArray after the last non-0 value
    histArrayShort = histArray(1:maxNonZeroLinks);
    
    hubindex = find(kpercent > 60) %Claim that hubs are any cell with >85% of connections
    nonhubindex = find(kpercent <= 60)
    highestconn = 0;
    for k = 1:length(hubindex)
        highestconn = [highestconn, find(N == hubindex(k))'];
    end
    lowconn = 0;
    for k = 1:length(nonhubindex)
        lowconn = [lowconn, find(N == nonhubindex(k))'];
    end
    highestconn(1) = [];
    lowconn(1)=[];
    [cellnum] = length(Rvars(i).data(:,10,l));

    
    Hubs(i).seed(l).data = highestconn;
    Hubsindex(i).seed(l).data = hubindex;
    Hubssum(l,i) = length(highestconn);
    linknum(l,i) = mean(hubindex);
   
    Kglychub(l,i) = (mean(Rvars(i).data(highestconn,10,l))-mean(Rvars(i).data(:,10,l)))/mean(Rvars(i).data(:,10,l));
    disp =  ((Rvars(i).data(highestconn,10,l))-mean(Rvars(i).data(:,10,l)))/mean(Rvars(i).data(:,10,l))
    
    end
end

%% Plot number of hubs cells
if Opts.figs
Hubcheck = 1;
fig = fig+1
figure(fig)
[h] = multiseedplot6(fig,Hubssum')
set(gca,'XTickLabel',Names);
title({'Number of Hubs (Hub has > 60% of cell connection)'})
ylabel({'Number of Hubs'})


fig = fig+1
figure(fig)

[h] = multiseedplot6(fig,Kglychub')
set(gca,'XTickLabel',Names);
title({'Percent difference Klgyc in hubs vs Islet'})
ylabel({'[Kglyc(hubs)-Kglyc(Islet)]/Kglyc(Islet)'})

end
%% 
if Opts.figs
if Hubcheck 
for l = 1:length(dataset2)
fig = fig +1;
[h] = multiseedplotsave(fig,squeeze(Pr_syn_gvn_gjhub(l,:,:))', ...
        ['Probability of Synchronization given Gap Junction' Names{l}], ...
        {'Hubs (>60% synchronization)', 'Low Links (< 60% synchronization)'},...
        {'Pr(Sync|GJ)'}, [savename '_Prsyncgvngjhubs' Names{l}])

fig = fig +1;
[h] = multiseedplotsave(fig,squeeze(Pr_gj_gvn_synhub(l,:,:))', ...
        ['Probability of Gap Junction given Synchronization' Names{l}], ...
        {'Hubs (>60% synchronization)','Low Links (< 60% synchronization)'},...
        {'Pr(GJ|Sync)'}, [savename '_Prgjgvnsynchubs' Names{l}])
    

fig = fig +1;
[h] = multiseedplotsave(fig,[Pr_syn_gvn_gj(l,:); Pr_gj_gvn_syn(l,:)], ...
        ['Probabilities ' Names{l}], ...
        {'Pr(Sync|GJ)', 'Pr(GJ|Sync)'},...
        {'Probability'}, [savename '_ProabilityAll' Names{l}])

end
end
end
%% 
if freqcalc == 1
for l = 1:5
       ca_freq = ca(1).data(2000:end, :,l);
       max_ca = max(ca_freq);
       bina = ca_freq./max_ca;
       bina_ac  = bina > .70;
       bina_fre = bina > .6;
       bina_in = bina < .4;       
       changeca = logical(diff(bina_fre));  

       for i = 1:cellnum
           freq(i) = nnz(changeca(:,i));
       end
       freq = mean(freq./2);

       %freq1000 = 1000*freq./size(ca_freq,1);
       period = size(ca_freq,1)./freq;
       periodsec(l) = period*100/1000;
end
if Opts.figs
fig = fig+1;
figure(fig), bar(periodsec)
title('Period of Oscillations')
xlabel('Seeds')
ylabel('Period (s)')
end
end
%% 
% fig = fig+1
% figure(fig)
%  [h] = multiseedplotsave(fig,Pr_both2,tit,xlab,ylab, savename)
% [h] = multiseedplot6(fig,Pr_both2)
% set(gca,'XTickLabel',dataset2);
% title({'Probabilty of functional connection given GJ - Second Method'})
% xlab = dataset2;
%xlab = string(Threshold)
xlab = Names
if Opts.figs
fig = fig+1
figure(fig)
% [h] = multiseedplot6(fig,Pr_gj_gvn_syn)
% set(gca,'XTickLabel',dataset2);
% title({'Probabilty of GJ given sync'})
if Hubcheck
tit = ['Probabilty of GJ given sync']
[h] = multiseedplotsave(fig,Pr_gj_gvn_syn,tit,xlab,'Pr(GJ|Sync)', ['PRGJgvnsync'])
end
fig = fig +1 ;
figure(fig)
if Hubcheck

tit = ['Probability of both GJ and Sync']
[h] = multiseedplotsave(fig,Pr_both, tit,xlab,'Pr(both)', ['PRboth' Names{l}])
end
fig = fig +1 ;
figure(fig)
if Hubcheck 

tit = ['Probability of both GJ']
[h] = multiseedplotsave(fig,Pr_gj, tit,xlab,'Pr(gj)', ['PRGJ' Names{l}])
end
fig = fig +1 ;
figure(fig)
if Hubcheck 

tit = ['Probability of Sync']
[h] = multiseedplotsave(fig,Pr_syn, tit,xlab,'Pr(sync)', 'PRsync')

end

fig = fig+1
figure(fig)
[h] = multiseedplot6(fig,Pr_syn_gvn_gj)
set(gca,'XTickLabel',Names);
title({'Probabilty of functional connection given GJ'})
% tit=({'Probabilty of functional connection given GJ'})
% [h] = multiseedplotsave(fig,Pr_syn_gvn_gj,tit,xlab,'Pr(Sync|GJ)', 'PRsyncgvnGJ')

fig = fig+1
figure(fig)
% 
[h] = multiseedplot6(fig,Net.Cavg)
set(gca,'XTickLabel',Names);
title({'Clustering Coeffiecient'})
% tit = 'Clustering Coeffiecient'
% [h] = multiseedplotsave(fig,Net.Cavg,tit,xlab,'Clustering Coefficient', 'CCoef')
% 


fig = fig+1
figure(fig)
[h] = multiseedplot6(fig, Net.L)
set(gca,'XTickLabel',Names);
title({'Caracteristic Path Length (excluding pairs with no path)'})
% tit = 'Caracteristic Path Length (excluding pairs with no path)';
% [h] = multiseedplotsave(fig,Net.L,tit,xlab,'Path Length', 'PathLength')


fig = fig+1
figure(fig)
[h] = multiseedplot6(fig,Net.Eglob)
set(gca,'XTickLabel',Names);
title({'Global Efficiency'})
% tit = 'Global Efficiency';
% [h] = multiseedplotsave(fig,Net.Eglob,tit,xlab,'Efficiency', 'Efficiency')


fig = fig+1
figure(fig)
[h] = multiseedplot6(fig, S)
set(gca,'XTickLabel',Names);
title({'Small-world-ness'})
% tit = 'Small-world-ness';
% [h] = multiseedplotsave(fig,S,tit,xlab,'Small World-ness', 'Small World-ness')



fig = fig+1
figure(fig)
[h] = multiseedplot6(fig, Lsm)
set(gca,'XTickLabel',Names);
title({'Path length/Random network path length'})
% tit = 'Path length/Random network path length'
% [h] = multiseedplotsave(fig,Lsm,tit,xlab,'L/L_{rand}', 'PathLengthRandom')
% 


fig = fig+1
figure(fig)
[h] = multiseedplot6(fig, Csm)
set(gca,'XTickLabel',Names);
title({'Clustering Coefficient/Random network Clustering Coefficient'})
% tit = 'Clustering Coefficient/Random network Clustering Coefficient'
% [h] = multiseedplotsave(fig,Csm,tit,xlab,'C/C_{rand}', 'CCRandom')


fig = fig+1
figure(fig)
[h] = multiseedplot6(fig, Esm)
set(gca,'XTickLabel',Names);
title({'Efficiency/Random network efficiency'})
% %
% tit = {'Efficiency/Random network efficiency'}
% [h] = multiseedplotsave(fig,Esm,tit,xlab,'E/E_{rand}', 'ERandom')

fig = fig+1
figure(fig)
[h] = multiseedplot6(fig, Csm)
set(gca,'XTickLabel',Names);
title({'Efficiency/Random network efficiency'})

fig = fig +1
figure(fig)

[h] = multiseedplot6(fig,Net.nopath)
set(gca,'XTickLabel',Names);
title({'% of pairs with no path'})
% tit = '% of pairs with no path';
% [h] = multiseedplotsave(fig,Net.nopath,tit,xlab,'% pairs with no path', 'no path')
end
%% plot phase
    fig = fig+1
% try
% phaseinfo.rand(3)
%  Isizes2 = Isizes
% catch
% Isizes2 = 2;
% end
if 0
for l = 1:Isizes
    for i = 1:seeds
    [cellnum] = length(Rvars(l).data(:,10,i));
    numcell = cellperc*cellnum;
    randomcell = randi(cellnum,numcell,1);
    N = Links2.N(l).data(:,i);

    phaseinfo.high(l).data(:,i) = mean(N(phaseinfo.cell(l).data(end-numcell+1:end,i)));
    phaseinfo.rand(l).data(:,i) = mean(N(randomcell));
    phaseinfo.low(l).data(:,i) = mean(N(phaseinfo.cell(l).data(1:numcell,i)));
    
    end
    if Opts.figs
    figure(fig)
    subplot(Opts.subplotdim+l)
    [h] = multiseedplot6(fig,[phaseinfo.high(l).data; phaseinfo.rand(l).data; phaseinfo.low(l).data])
    set(gca,'XTickLabel',{'High Phase (lagging)', 'Rand Links', 'Low Phase (leading)'});
    title({'Average number of links' '(Islet Size = ' num2str(cellnum) 'cells)' Names{l}})
    ylabel({'Links'})
    end
end
end



save('/Volumes/Briggs_2TB/SizeDependence/Analysis/Coup_Newclustering_savetime.mat')
