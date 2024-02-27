%This is the main script which loads all of the calcium data and runs
%analyses used in Briggs et al., 2023
% Jennifer Briggs, 2021

clear all
close all



%%%%%%%%%%%%%%%%%%%%%%%%%%
%Change current working directory to the Functional_and_Structural_Networks
%Parent folder
cd = pwd;
addpath(genpath(cd)); %adding all subfolders to path

savetime =  datestr(datetime('today'),'yyyymmdd');

%%This was just the file structure that I used to work with different
%%experiments and trials
dataset2 = {'ExampleData/simulation/FastOscillations/'}
dataname = {'/seed1/', '/seed2/'}; %these are your repeat trials
datapath = pwd


cellperc = .1; %percentage of cells to look at (if you are looking at top x% or something)

%You can either set the threshold or automatically choose a threshold using
%optRth.m
Manual_Th = 0; %set to 1 if you want to choose a threshold, 0 will call optRth.m
Threshold = .9995; %If Manual_Th = 1, you need to set a threshold. Otherwise this will be overwritten

savename = ['savename' savetime]; 

gjconnectionprob = 1 %define whether you want to compare the functional network to gap junctions (1) or metabolism (0)

%Set options for network analysis
Opts.figs = 0;
Opts.change_thres = 0;
Opts.chage_phase = 0;
Opts.printasSubplot = 1;
Opts.subplotdim = 310; %Set the last number to 0 always, first number is rows second is columns
Opts.multiseed = 1;

fig = 1; %counter for figures

%define data lengths
seeds = length(dataname);
Isizes = length(dataset2);
%% Import Data

%make data structures for calcium (ca) and parameters (such as katp, kglyc,
%etc, (Rvars). 
ca = struct; 
Rvars = struct;

%% Loading Data %%
tic
for l = 1:length(dataset2)
    for i = 1:length(dataname)
    disp(['Importing ' cell2mat(dataset2(l))])
    capath = fullfile(datapath, dataset2(l), dataname(i), 'calcium.txt'); 
    casize = importdata(cell2mat(capath));

    %this is just a big loop to help ensure that all of the calcium arrays
    %are equal sizes
    if i == 1
       ca2(l).data(:,:,i) = casize;
    else
       sizes = [size(casize,1), size(ca2(l).data(:,:,i-1),1)]
       [sizemin,seedmin] = min(sizes);
       for tt = 1:l
       ca2(tt).data(sizemin+1:end,:,:)=[]
       end
       casize(sizemin+1:end,:)=[];
       ca2(l).data(:,:,i) = casize;
    end

    %extracting information
    ca2(l).name = dataset2(l);
    
    %positions
    pospath = fullfile(datapath, dataset2(l), dataname(i), 'XYZpos.txt');
    pos(l).data(:,:,i) = importdata(cell2mat(pospath));
    pos(l).name = dataset2(l);

    %Gap Junction Connections
    connpath = fullfile(datapath, dataset2(l), dataname(i), 'NN10A.txt');
    %this variable contains the gap junction connection structure in the
    %simulation. Each row (i) corresponds to cell (i) and each column contains
    %the index of the other cells that cell (i) is connected to. -1 is the
    %filler for cells without connections. BUT this is written in C, so the
    %indexing in matlab is the index +1
    conn(l).data(:,:,i) = importdata(cell2mat(connpath));  

    % extract second phase of calcium
    calcium = ca2(l).data(:,:,i);    
    ca(l).data(:,:,i) = calcium(2000:end,:); %second phase only

    % Extract parameters such as Kglyc, Katp, etc. 
    varspath = fullfile(datapath, dataset2(l), dataname(i), 'RandomVars.txt');
    Rvars(l).data(:,:,i) = importdata(cell2mat(varspath));
    Rvars(l).name = dataset2(l);
    fig1 = figure(1)%, subplot(Opts.subplotdim+l)
    
    
    nexttile
    plot(ca(l).data(:,:,i)), xlabel('Time'), ylabel('Ca'), title(dataname(i))
    
    
   
    end
end
toc

clear casize sizes sizemin seed min tt


%% Run Network Analysis
Links2 = struct;
Links2.N = struct;
Links2.sort = struct;
Links2.cellsort = struct;
fig = fig+1;
for l = 1:length(dataset2)
    for i = 1:length(dataname)
        
     Opts.printasSubplot = 1;  Opts.Subplotnum = l; Opts.multiseednum = i; %define options
         
     if Manual_Th == 0 %if the threshold is not set, we automatically define one
         Opts.avDeg = 3
         Opts.Min = 0
         Opts.Method = 'Degree'
         Threshold = double(findoptRth(ca(l).data(:,:,i), Opts));
     end

     %run links analysis
    [Links2.N(l).data(:,i), adjm(l).data(:,:,i), kperc(l).Islet(i).data , histArrayPercShort(i).data, Rij(l).Islet(i).data(:,:) ] = links(ca(l).data(:,:,i), Threshold,Opts,fig); %%This is where network analysis is perfromed

        %sort
       [Links2.sort(l).data(:,i), Links2.cellsort(l).data(:,i)]= sort(Links2.N(l).data(:,i));
    end
end


%% Here is where we look at the parameter values corresponding to the functional network (mostly figure 1)
%predefine objects
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
    %Sort cells in ascending order based on Kglyc. 'sort' is the parameter
    %value and 'cell' is the corresponding cell index.
    [Kglyc.sort(l).data(:,i), Kglyc.cell(l).data(:,i)]=sort(Rvars(l).data(:,10,i));
    Kglyc.kglyc(l).data(:,i) = (Rvars(l).data(:,10,i)); %unsorted parameter value



    [cellnum] = length(Rvars(l).data(:,10,i));%just extracting how many cells there are.
    numcell = cellperc*cellnum; %cellperc was set above for looking at the top x% of cells. 
    randomcell = randi(cellnum,numcell,1); %select random 10% of cells

    %extract the degree of each cell
    N = Links2.N(l).data(:,i);
    cell = Kglyc.cell(l).data(:,i);


    %extract the functional degree for the high, low, and random cells
    %sorted by their kglyc value
    Kglyc.high(l).data = mean(N(cell(end-numcell:end)));
    Kglyc.rand(l).data = mean(N(cell(randomcell)));
    Kglyc.low(l).data = mean(N(cell(1:numcell)));

    %actual degree highest, random, and lowest degree cells. Note that
    %since we are working with a scale-free-like distribution, there are a
    %lot of 0 degree cells. 
    Links2.high(l).data(:,i) = (Links2.cellsort(l).data(end-numcell:end,i));
    Links2.rand(l).data(:,i) = (Links2.cellsort(l).data(randomcell,i));
    Links2.low(l).data(:,i) = (Links2.cellsort(l).data(1:numcell,i));
    
    % We just extracted cells/param values based on percentages but we can also
    % Find hubs as defined by Johnston et al. 2016 
    a = find(kperc(l).Islet(i).data  > 60);
    hubthreshold = (a(1));
    b = find(kperc(l).Islet(i).data  <= 60);
    connect = sum(adjm(l).data(:,:,i)); 
    Hubs(l).hubby(i).data = find(connect>hubthreshold);%Number of cells linked to more than 60% of Islet
    hubby = Hubs(l).hubby(i).data;
    
   %looking at non hub cells
	hublow = find(connect<=hubthreshold);
    Hubs(l).hublow(i).data = hublow;%Number of cells linked to more than 50% of Islet

    %extract kglyc values for hubs/con hubs
    Kglyc.hubs(l).data(:,i) = mean(Kglyc.kglyc(l).data(hubby,i));
    Kglyc.hubslow(l).data(:,i) = mean(Kglyc.kglyc(l).data(hublow,i));

    %extract kglyc values based on top x% of cells rather than hubs
    Kglyc.high(l).data(:,i) = mean(Kglyc.kglyc(l).data(Links2.high(l).data,i));
    Kglyc.rand(l).data(:,i) = mean(Kglyc.kglyc(l).data(Links2.rand(l).data,i));
    Kglyc.low(l).data(:,i) = mean(Kglyc.kglyc(l).data(Links2.low(l).data,i));
   

    %extract STRUCTURAL network
    Conny = sum((conn.data(:,:,i)~=-1)');
    
    TotConn.hubs(l).data(:,i) =mean(Conny(hubby)); %average number of physical connections shared by hubs
    TotConn.hubslow(l).data(:,i) = mean(Conny(hublow)); %average number of physical connections shared by non hubs
    
    Katp.hubs(l).data(:,i) = mean(conn.data(hubby,1,i));
    Katp.hubslow(l).data(:,i) = mean(conn.data(hublow,1,i));

    
    % average GJ conductance is the sum of average gj conductance between
    % connected cells
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
     
    %Find the gap junction conductance for hubs/non hubs
    Gcoup.hubs(l).data(:,i) = mean(GJconduct(i).data(hubby));
    Gcoup.hubslow(l).data(:,i) = mean(GJconduct(i).data(hublow));

    %extract kglyc values based on top x% of cells rather than hubs
    Gcoup.high(l).data(:,i) = mean(Rvars(l).data(Links2.high(l).data,2,i));
    Gcoup.rand(l).data(:,i) = mean(Rvars(l).data(Links2.rand(l).data,2,i));
    Gcoup.low(l).data(:,i) = mean(Rvars(l).data(Links2.low(l).data,2,i));
end

    
   clear N cell
end

 for i = 1:size(dataname)
    [~,sorted] = sort(Links2.N.data(:,i));
    Gcoup.top10(i) = mean(GJconduct(i).data(sorted(end-numcell:end)))
    Gcoup.bottom90(i) = mean(GJconduct(i).data(sorted(1:end-numcell)))
    kglycdata = Rvars.data(:,10,i);
    katpdata = Rvars.data(:,1,i);
    Kglyc.top10(i) = mean(kglycdata(sorted(end-numcell:end)));
    Kglyc.bottom90(i) = mean(kglycdata(sorted(1:end-numcell)));
    Katp.top10(i) = mean(katpdata(sorted(end-numcell:end)));
    Katp.bottom90(i) = mean(katpdata(sorted(1:end-numcell)));
end

%% FIGURE 4: Probability that a cell is linked given functional connection
clear t counter
for i = 1:Isizes
   for l = 1:seeds
       tic

    posit = pos(i).data(:,:,l); %XYZ data

    %get Structural connections
    gjconn = conn(i).data(:,:,l);
    gjconn = gjconn + 1; % Data is in C++ so we need to start with 1 instead of 0

    %average number of structural connections
    avgjconn(i,l) = sum(nnz(gjconn))/size(gjconn,1);

    %functional connections
    adj = adjm(i).data(:,:,l);
    adj = adj - diag(diag(adj)); % Make adj matrix diags 0

  %Convert gjconn to an adjacency matrix (structural network)
   gj2adj = zeros(size(adj));
        t = 1;
        if gjconnectionprob == 1 %Create a network of gap junction connected cells
            for c = 1:size(gjconn,1)
                for n=1:size(gjconn,2)
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

    if 1
    dist = sqrt((posit(:,1)-posit(:,1)').^2+(posit(:,2)-posit(:,2)').^2+...
    (posit(:,3)-posit(:,3)').^2);
    dist(dist==0)=100;
    end

    %threshold for similar kglyc
    target = mean2(Rvars(i).data(:,10,l)) + 1*std(Rvars(i).data(:,10,l)); 

if gjconnectionprob == 0 %Create a network of kglyc cells > some value
    disty = []

    for c = 1:cellnum %for each cell, 
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
         for n=1:cellnum
                gj2adj(c,n) = abs((Rvars(i).data(c,10,l)-Rvars(i).data(n,10,l)));     %setting up the matrix of the difference between each kglyc  
         end  
    end
 
  gj2adj(dist > .15) = 0; %distance has to be restricted.

  %draw a connection based on if kglyc is more similar than islet average
  gj2adj = gj2adj<(mean2(nonzeros(gj2adj))-0*std(nonzeros(gj2adj),0,'all')) & gj2adj>0;

  gj2adj = double(gj2adj); %now the adjacency graph is binary
end
% 

%% FIGURE 5: shortest path of the weighted graph :
    mm(l) = mean(sum(gj2adj, 'omitnan'));
    
    if gjconnectionprob
        gjconduct = Rvars(i).data(:,10,l);  %Weight by GJ
    else
        gjconduct = Rvars(i).data(:,2,l); %Weight by Metabolism
    end
    [distweighted(:,:,l), adjw] = weightednet(gj2adj, gjconduct,0); %outputs graph weighted by inverse of conductance
    [Length, syncD_all, nonsyncD_all] = WeightedDistEqualbasedonCell(adj,adjw,string(l)); %calculates the weight and distance for syncrhonized and non syncrhonized cells


    %If you don't care about separating by distances, just calculate the
    %weighted distances for all syncrhonzied cells
    [col,row] = find(adj); %finding edges
    dd = distweighted(:,:,l);
    dd(isinf(dd))= NaN; %NaN if there is no path between two cell pairs
    dd(dd==0) = NaN;
    for jj = 1:length(col)
        ShortestpathSync(jj) = dd(col(jj),row(jj));
    end
    numcells = length(adj);
     
    x=[];
    y=[];


%% experimental part of Figure 7: 
        [Net.L(i,l),Net.Lrand(i,l), Net.Eglob(i,l), ... 
         Net.Erand(i,l), Net.Eloc(i,l), Net.nopath(i,l),...
         ~, Net.Cavg(i,l), Net.Crand(i,l)] = graphProp(adj,x,y,0)
        S_stz(i,l) = (Net.Cavg(i,l)/Net.Crand(i,l))/(Net.Erand(i,l)/Net.Eglob(i,l));
        Lsm(i,l) =  Net.L(i,l)/Net.Lrand(i,l);
        Csm(i,l) = Net.Cavg(i,l)/Net.Crand(i,l);
        Esm(i,l) = Net.Erand(i,l)/Net.Eglob(i,l);
        S_strog(i,l) = (Net.Cavg(i,l)/Net.Crand(i,l))/(Net.L(i,l)/Net.Lrand(i,l));



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


    Pr_both(i,l) = (sum(nnz(both))/2)/(((size(both,1)-1)*(size(both,1))/2));

    %if gjconnectionprob == 0 then this is actual probability of
    %syncrhonization given *metabolic* similarity and probability of
    %*metabolic* similarity given synchronization. Otherwise it is gap
    %junction and synchronization
    Pr_syn_gvn_gj(i,l) = Pr_both(i,l)/Pr_gj(i,l);
    Pr_gj_gvn_syn(i,l) = Pr_both(i,l)/Pr_syn(i,l); 

    Cpos =   (((size(adj,1)-1)*(size(adj,1))/2));
    Links2.Ci(i).data(:,l) = Links2.N(i).data(:,l)./Cpos;

    clear counter
    toc
    end
end

