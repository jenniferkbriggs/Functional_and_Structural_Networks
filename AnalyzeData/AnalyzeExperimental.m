%% Analyze experimental data:
% This code is written to analyze networks for experimental data from Anne
% and Vira. Files are calcium data.

%In this script, we load and analyze all of the imaging files for a single analysis 
addpath('AnalyzeData/')

filesdir = ('ExampleData/experimental/'); %directory of your files
files = dir([filesdir '*CA.csv']); %list all files in directory
ct = 1;

%%first find optimal threshold defined as the average optimal threshold for
%%all imaging files
for i = 1:length(files) %indicies correspond to file indices to sort through
    CA = readmatrix([filesdir files(i).name]);
    time = CA(:,1);
    CA(:,1) = [];
    figure, plot(CA), title('Identify beginning of second phase') % cut out first phase if applicable.
    [st(ct), ~] = ginput(1)
    CA = CA(st(ct):end, :);
    close all


    %here we chose what method to use to identify optimal threshold. See documentation in "findoptRth"
    Opts.Method = 'Scale-Free'; 
    Opts.Max = 10;
    Opts.Min = 3;
    thr(ct) = findoptRth(CA, Opts)
     ct = ct+1
end

%% now find hubs and assign degree:
ct = 1
for i = 1:length(files)
    CA = readmatrix([filesdir files(i).name]);
    CA = CA(st(ct):end, :);

    % if size(CA,1)<size(CA,2)
    %     CA = CA';
    % end

    Rij = corr(CA);
    Adj = (Rij > thr);
    degree(1:length(Rij),ct) = sum(Adj);

    normdegree(1:length(Rij),ct) = sum(Adj)./max(sum(Adj));
    numhubs(ct) = length(find(normdegree(:,ct)>= 0.6))./length(Adj);
    ct = ct+1
end

hubs = nan(size(degree));
for i = 1:length(files)
    normdegree(:,i) = degree(:,i)./max(degree(:,i));
    numhubs(i) = length(find(normdegree(:,i)>= 0.6))./length(degree(:,i));
    hubs(1:length(find(normdegree(:,i)>= 0.6)),i) = find(normdegree(:,i)>= 0.6); %find cell index of hubs
end


%% Sort FRAP by degree (Figure 3g)
% - every dot indicates all FRAP greater than a certain degree and less than the next degree. So we
% don't double count cells
FRAP = readmatrix([filesdir 'exampleFRAP.csv'])
for i = 1:size(FRAP,2) %if using example data, there is only one islet
    indx = find(~isnan(FRAP(:,i))); %find the cells that we also have FRAP data for
    degrees = normdegree(indx,i);
    frap = FRAP(indx,i); %get the frap data
    degrees = floor(degrees*10); %convert the degrees into which bin they will be located in

    [sorteddeg, sortindx]= sort(degrees)

    %change 10 to 9 because we bin by > 0.9
    if length(find(sorteddeg == 10)) > 0
        sorteddeg(find(sorteddeg == 10)) = 9
    end

    %sort frap by degree
    frap_sorted = frap(sortindx);
    udeg = unique(sorteddeg);

    
    %count how many cells are in each bin
    for j = 1:length(udeg)
        ct = (udeg(j))    
        %find frap indices that correspond to unique degree
        allindx = (find(sorteddeg == ct));

        frap_final(i,ct+1) = mean(frap_sorted(allindx));
    end
end
degreebins = [9:-1:0]; %>0.9	>0.8	>0.7	>0.6	>0.5	>0.4	>0.3	>0.2	>0.1	<0.1
frap_final = fliplr(frap_final); %this is just because we plotted the histogram reversed from how it is calculated (ascending rater than descencing)


%% Get just the hub cells and top10% of cells = 
% - every dot indicates all FRAP greater than a certain degree and less than the next degree. So we
% don't double count cells
FRAP = readmatrix([filesdir 'exampleFRAP.csv'])
hubs_frap = nan(22,size(FRAP, 2));
nonhubs_frap = nan(33,size(FRAP, 2));

%50 is the maximum number of photobleached cells in our dataset (not
%necessarily the example csv)- this needs to be altered to fit your dataset
top10_frap = nan(50,size(FRAP, 2));
bottom90_frap = nan(50,size(FRAP, 2));
for i = 1:size(FRAP, 2)

    %sort top 10% of cells and bottom 90%
        cutoff = sort(normdegree(:,i));
        cutoff = cutoff(~isnan(cutoff));
        cutoff = cutoff(end-floor(0.1*length(cutoff)):end);
        top10indx = find(normdegree(:,i) >= cutoff(1));
        bottom90indx = setxor([1:length(find(~isnan(degree(:,i))))], top10indx);

        top10_frap(1:length(top10indx),i) = (FRAP(top10indx,i));
        bottom90_frap(1:length(bottom90indx),i) = (FRAP(bottom90indx,i));


    %find hubs (<60% of max degree)
        hubsindx = hubs(~isnan(hubs(:,i)),i);
        hubs_frap(1:length(hubsindx),i) = (FRAP(hubsindx,i));
        
        nhubindx = setxor([1:length(find(~isnan(degree(:,i))))], hubsindx);
        nonhubs_frap(1:length(nhubindx),i) = (FRAP(nhubindx,i));
        
end

nonhubsfrap = nonhubs_frap(~isnan(nonhubs_frap));
hubsfrap = hubs_frap(~isnan(hubs_frap));

top10frap = top10_frap(~isnan(top10_frap));
bottom90frap = bottom90_frap(~isnan(bottom90_frap));



