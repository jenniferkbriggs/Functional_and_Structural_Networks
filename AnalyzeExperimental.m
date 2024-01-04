%% Analyze experimental data:
% This code is written to analyze networks for experimental data from Anne
% and Vira. Files are calcium data.


filedir = '/Volumes/Briggs_10TB/NetworkPaper/Experimental_Data/Slow_oscillations/'
data = dir(filedir)
ct = 1;

% %first find optimal threshold
for i = 12:20 %indicies correspond to file indices to sort through

    CA = readmatrix([filedir data(i).name]);
    figure, plot(CA), title('Identify beginning of first response')
   [st(ct), ~] = ginput(1)
    CA = CA(st(ct):end, :);
    close all

    Opts.Method = 'Scale-Free';
    Opts.Max = 10;
    Opts.Min = 3;
    thr(ct) = findoptRth(CA, Opts)
     ct = ct+1
end

thr = 0.9
ct = 1
%% now find hubs and assign degree:
for i = 12:20
    CA = readmatrix([filedir data(i).name]);
    CA = CA(st(ct):end, :);

    % if size(CA,1)<size(CA,2)
    %     CA = CA';
    % end

    Rij = corr(CA);
    Adj = (Rij > thr);
    degree(1:length(Rij),ct) = sum(Adj);

    normdegree(1:length(Rij),ct) = sum(Adj)./max(sum(Adj));
    numhubs(ct) = length(find(normdegree(:,ct)>= 0.6))./length(Adj);
    degree(length(Rij)+1:50,ct) = NaN
    ct = ct+1
end

hubs = nan(size(degree));
for i = 1:9
    normdegree(:,i) = degree(:,i)./max(degree(:,i));
    numhubs(i) = length(find(normdegree(:,i)>= 0.6))./length(degree(:,i));
    hubs(1:length(find(normdegree(:,i)>= 0.6)),i) = find(normdegree(:,i)>= 0.6);
end


%% My method
% - every dot indicates all FRAP greater than a certain degree and less than the next degree. So we
% don't double count cells
FRAP = readmatrix(['/Volumes/Briggs_10TB/NetworkPaper/Experimental_Data/FRAP.xlsx'])
frap_final = nan(10, 10);
for i = 1:9
    indx = find(~isnan(FRAP(:,i)));
    degrees = normdegree(indx,i);
    frap = FRAP(indx,i)
    degrees = floor(degrees*10);
    [sorteddeg, sortindx]= sort(degrees)

    %change 10 to 9 because we bin by > 0.9
    if length(find(sorteddeg == 10)) > 0
        sorteddeg(find(sorteddeg == 10)) = 9
    end
    frap_sorted = frap(sortindx);
    udeg = unique(sorteddeg);

    
    
    for j = 1:length(udeg)
        ct = (udeg(j))    
        %how many?
        allindx = (find(sorteddeg == ct));

        frap_final(i,ct+1) = mean(frap_sorted(allindx));
    end


end

frap_final = fliplr(frap_final);


%% Get just the hub cells = 
% - every dot indicates all FRAP greater than a certain degree and less than the next degree. So we
% don't double count cells
FRAP = readmatrix(['/Volumes/Briggs_10TB/NetworkPaper/Experimental_Data/FRAP.xlsx'])
FRAP(50,:) = NaN;
hubs_frap = nan(22,9);
nonhubs_frap = nan(33,9);

top10_frap = nan(50,9);
bottom90_frap = nan(50,9);
for i = 1:9
        cutoff = sort(normdegree(:,i));
        cutoff = cutoff(~isnan(cutoff));
        cutoff = cutoff(end-floor(0.1*length(cutoff)):end);
        top10indx = find(normdegree(:,i) >= cutoff(1));
        bottom90indx = setxor([1:length(find(~isnan(degree(:,i))))], top10indx);

        top10_frap(1:length(top10indx),i) = (FRAP(top10indx,i));
        bottom90_frap(1:length(bottom90indx),i) = (FRAP(bottom90indx,i));


%find hubs
        hubsindx = hubs(~isnan(hubs(:,i)),i);
        hubs_frap(1:length(hubsindx),i) = (FRAP(hubsindx,i));
        
        nhubindx = setxor([1:length(find(~isnan(degree(:,i))))], hubsindx);
        nonhubs_frap(1:length(nhubindx),i) = (FRAP(nhubindx,i));
        
end

nonhubsfrap = nonhubs_frap(~isnan(nonhubs_frap));
hubsfrap = hubs_frap(~isnan(hubs_frap));

top10frap = top10_frap(~isnan(top10_frap));
bottom90frap = bottom90_frap(~isnan(bottom90_frap));
%% Vira's method
% - every dot indicates all FRAP greater than a certain degree. So we
% double count cells. THIS IS USED IN THE FINAL MANUSCRIPT AS IT WAS THE
% METHOD PRESENTED IN THE ORIGINAL SUBMISSION
FRAP = readmatrix(['/Volumes/Briggs_10TB/NetworkPaper/Experimental_Data/FRAP.xlsx'])
frap_final = nan(10, 10);
for i = 1:9
    indx = find(~isnan(FRAP(:,i)));
    degrees = normdegree(indx,i);
    frap = FRAP(indx,i)

    for j = 9:-1:6
       allindx = (find(degrees >= j/10)); 
       frap_final(i,j+1) = mean(frap(allindx));
    end

    for j = 6:-1:1
    allindx=find(degrees<=j/10);
    frap_final(i,j) = mean(frap(allindx));
    end


end  


