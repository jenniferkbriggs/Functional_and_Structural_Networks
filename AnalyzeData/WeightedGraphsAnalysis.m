%This script was used to analyze the networks based on weight - given in
%figure 5.

type = ''
one = load([type 'ShortPathEqualDist.mat']);
two = load([type 'ShortPathEqualDist2.mat']);
three = load([type 'ShortPathEqualDist3.mat']);
four = load([type 'ShortPathEqualDist4.mat']);
five = load([type 'ShortPathEqualDist5.mat']);

%get one islet:
nonSyncone = NaN([20,300000]);
Syncone = NaN([20,300000]);

for i = 1:20
    nonSyncone(i,1:length(find(one.Length==i))) = 1./(one.nonsyncD_all(find(one.Length == i)));
    Syncone(i,1:length(find(one.Length==i)))  = 1./(one.syncD_all(find(one.Length == i)));
end

for i = 1:20
    nonSync(i,1) = 1./mean(one.nonsyncD_all(find(one.Length == i)),'omitnan');
    Sync(i,1)  = 1./mean(one.syncD_all(find(one.Length == i)),'omitnan');
    nonSync(i,2) = 1./mean(two.nonsyncD_all(find(two.Length == i)),'omitnan');
    Sync(i,2)  = 1./mean(two.syncD_all(find(two.Length == i)),'omitnan');
    nonSync(i,3) = 1./mean(three.nonsyncD_all(find(three.Length == i)),'omitnan');
    Sync(i,3)  = 1./mean(three.syncD_all(find(three.Length == i)),'omitnan');
    nonSync(i,4) = 1./mean(four.nonsyncD_all(find(four.Length == i)),'omitnan');
    Sync(i,4)  = 1./mean(four.syncD_all(find(four.Length == i)),'omitnan');
    nonSync(i,5) = 1./mean(five.nonsyncD_all(find(five.Length == i)),'omitnan');
    Sync(i,5)  = 1./mean(five.syncD_all(find(five.Length == i)),'omitnan');
end


for i = 1:20
    ynonSync(i,1) = mean(1./one.nonsyncD_all(find(one.Length == i)),'omitnan');
    ySync(i,1)  = mean(1./one.syncD_all(find(one.Length == i)),'omitnan');
    ynonSync(i,2) = mean(1./two.nonsyncD_all(find(two.Length == i)),'omitnan');
    ySync(i,2)  = mean(1./two.syncD_all(find(two.Length == i)),'omitnan');
    ynonSync(i,3) = mean(1./three.nonsyncD_all(find(three.Length == i)),'omitnan');
    ySync(i,3)  = mean(1./three.syncD_all(find(three.Length == i)),'omitnan');
    ynonSync(i,4) = mean(1./four.nonsyncD_all(find(four.Length == i)),'omitnan');
    ySync(i,4)  = mean(1./four.syncD_all(find(four.Length == i)),'omitnan');
    ynonSync(i,5) = mean(1./five.nonsyncD_all(find(five.Length == i)),'omitnan');
    ySync(i,5)  = mean(1./five.syncD_all(find(five.Length == i)),'omitnan');
end



for i = 1:20
    innonSync(i,1) = mean(one.nonsyncD_all(find(one.Length == i)),'omitnan');
    inSync(i,1)  = mean(one.syncD_all(find(one.Length == i)),'omitnan');
    innonSync(i,2) = mean(two.nonsyncD_all(find(two.Length == i)),'omitnan');
    inSync(i,2)  = mean(two.syncD_all(find(two.Length == i)),'omitnan');
    innonSync(i,3) = mean(three.nonsyncD_all(find(three.Length == i)),'omitnan');
    inSync(i,3)  = mean(three.syncD_all(find(three.Length == i)),'omitnan');
    innonSync(i,4) = mean(four.nonsyncD_all(find(four.Length == i)),'omitnan');
    inSync(i,4)  = mean(four.syncD_all(find(four.Length == i)),'omitnan');
    innonSync(i,5) = mean(five.nonsyncD_all(find(five.Length == i)),'omitnan');
    inSync(i,5)  = mean(five.syncD_all(find(five.Length == i)),'omitnan');
end