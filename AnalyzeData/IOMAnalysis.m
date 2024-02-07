%% Analyze IOM:
cd('/Volumes/Briggs_10TB/NetworkPaper/IOMSlow/')
%clear all 
close all
%predefine threshold with very high variable precision

if 0 %change to 0 to do variable precision
    stoc = 1;
    rng(3)
else
    %thr = vpa([]);
    stoc = 0
end

for i = 1:5
    Params = readmatrix(['IOM2_set11_250cell_' num2str(i) '_parameters.csv']); %cell, Kglyc, KATP
    Conduct = readmatrix(['IOM2_set11_250cell_' num2str(i) '_adjacency_matrix.csv']);
    Ca5 = readmatrix(['IOM2_set11_250cell_' num2str(i) '_calcium.csv']);
   
    time = Ca5(:,1);
    Ca5(:,1) = [];
    
    figure(1), nexttile
    plot(time, Ca5)
    title(['Set ' num2str(i)])

    
    %add a very small amount of noise to highlight the differences in
    %correlation

    if stoc
     noise = 0.001*rand(size(Ca5));
     Ca5 = noise + Ca5;
    end
    [Rij, pval]=corr((Ca5));
    Rij = (Rij)-(diag(diag(Rij)));
    Rij = vpa(Rij); %increase machine precision
    
    Opts.avDeg = 3
    Opts.Min = 0
    Opts.Method = 'Degree'
   %thr(i) = findoptRth(Ca5, Opts)

   %fast:
    %Adj =  double(Rij>0.9999993);
    %Adj = double(Rij>0.9999994)
   %slow:
    Adj =  double(Rij>0.999999999);
 %   Adj = double(Rij > thr(i));
   % Adj = double(Rij < 0.99999871513);

   [~,sorted] = sort(sum(Adj));

   numcell = length(Conduct);
    top10 = sorted(end-(numcell*.10):end);
    bottom90 = sorted(1:end-numcell*.1);

    Hubs = find(sum(Adj) > 0.6*max(sum(Adj)));
    Hub_indx(i,1:length(Hubs)) = Hubs;
    nonHubs = setxor([1:length(Adj)], Hubs);
    numhubs(i) = length(Hubs)

    %kglyc is 2nd
    hubs_kglyc(i,1:length(Hubs)) = Params(Hubs, 2);
    nonhubs_kglyc(i,1:length(nonHubs)) = Params(nonHubs, 2);
    hubs_katp(i,1:length(Hubs)) = Params(Hubs, 3);
    nonhubs_katp(i,1:length(nonHubs)) = Params(nonHubs, 3);
    
    Gcoup = sum(Conduct);
    hubs_Gcoup(i,1:length(Hubs)) = Gcoup(Hubs);
    nonhubs_Gcoup(i,1:length(nonHubs)) = Gcoup(nonHubs);

    degrees(i,:) = sum(Adj)./max(sum(Adj));

    top10_kglyc(i) = mean(Params(top10,2));
    bottom90_klgyc(i) = mean(Params(bottom90,2));
    top10_gcoup(i) = mean(Gcoup(top10));
    bottom90_gcoup(i) = mean(Gcoup(bottom90));
    top10_katp(i) = mean(Params(top10,3));
    bottom90_katp(i) = mean(Params(bottom90,3));
end

for i = 1:5
Nad_hub(1:length(find(~isnan(Hub_indx(:,i)))),i) = NAD(Hub_indx(~isnan(Hub_indx(:,i)),i));
end

for i = 1:5
kglyc(i,:) = [mean(nonzeros(hubs_kglyc(i,:))), mean(nonzeros(nonhubs_kglyc(i,:)))];
katp(i,:) = [mean(nonzeros(hubs_katp(i,:))), mean(nonzeros(nonhubs_katp(i,:)))];
gcoup(i,:) = [mean(nonzeros(hubs_Gcoup(i,:))), mean(nonzeros(nonhubs_Gcoup(i,:)))];
end

for i = 1:5
nonhub_indx(1:length(setxor([1:length(Adj)], Hub_indx(~isnan(Hub_indx(:,i))))),i) = setxor([1:length(Adj)], Hub_indx(~isnan(Hub_indx(:,i))));
end

for i = 1:5
nonhub_nad(1:length(setxor([1:length(Adj)], Hub_indx(~isnan(Hub_indx(:,i))))),i) = NAD(nonhub_indx(find(nonhub_indx(:,i)~=0),i),i)
end

% %Look at graph:
% G = graph(Conduct>0)
% figure, plot(G