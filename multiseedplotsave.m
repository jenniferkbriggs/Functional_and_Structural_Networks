function [h] = multiseedplotsave(n,data,tit,xlab,ylab, savename)
%Plots mutliseed bar chart with significance
%Jennifer Briggs 2020

%Inputs: 
%n = figure number
%data = 2d data array
%tit = title: string
%xlab = xlabel: string
%ylab = ylabel: string
%savename = name to save figure: string
savetime =  datestr(datetime('today'),'yyyymmdd');

figure(n)

datamean = [mean(data')];
b=bar(datamean, 'FaceColor','flat');
hold on

b.CData(1,:) = [1 1 0];
b.CData(2,:)= [0 1 0];
if length(datamean) == 3
b.CData(3,:) = [1 0 0];
elseif length(datamean) == 4
b.CData(3,:) = [1 0 0];
b.CData(4,:) = [1 0 1];
elseif length(datamean) == 5
b.CData(3,:) = [1 0 0];
b.CData(4,:) = [1 0 1];
b.CData(5,:) = [0 1 1];
elseif length(datamean) == 6
b.CData(3,:) = [1 0 0];
b.CData(4,:) = [1 0 1];
b.CData(5,:) = [0 1 1];
b.CData(6,:) = [1 .7 .7];
end

for i = 1:length(datamean)
    KG_table_size = [data(i,:)];
    plot(i,KG_table_size','*')
end
if i == 3
    
h1 = ttest2(data(1,:),data(2,:))
h2 = ttest2(data(1,:),data(3,:))
h3 = ttest2(data(2,:),data(3,:))

h = [h1 h2 h3];

ctr2 = bsxfun(@plus, b(1).XData, [b(1).XOffset]');
hold on
y = max(datamean);
if h1 == 1
    plot([1*1.1, 2*.9], [1 1]*y*1.03, '-k', 'LineWidth',2)
    %plot(mean([1*1.1, 2*.9]), y*1.05+.03, '*k')
end
if h2 ==1 
    plot(ctr2(1:3), [1 1 1]*y*1.07, '-k', 'LineWidth',2)
   % plot(mean(ctr2(1:3)), y*1.1+0.3, '*k')
end
if h3 ==1
    plot([2*1.1,3*.9], [1 1]*y*1.03, '-k', 'LineWidth',2)
    %plot(mean([2*1.1,3*.9]), y*1.05+0.3, '*k')
end

ylim([0.9*min(min(data)), 1.1*max(max(data))])
end
if i == 2
    h = ttest2(data(1,:),data(2,:))

ctr2 = bsxfun(@plus, b(1).XData, [b(1).XOffset]');
hold on
y = max(datamean);
if h == 1
    plot([1*1.1, 2*.9], [1 1]*y*1.03, '-k', 'LineWidth',2)
end
else 
    h = 0
end

ylim([0.9*min(min(data)), 1.1*max(max(data))])
set(gca,'XTickLabel',xlab);
title(tit)
ylabel(ylab)

set(gcf,'position',[10,20,1000,900])

saveas(gcf, [pwd '\' savename '_' savetime '.fig'])
saveas(gcf, [pwd '\' savename '_' savetime '.png'])

data = array2table(data', 'VariableNames',xlab)
writetable(data, [pwd '\' savename '_' savetime '.xlsx'])

save([pwd '\' savename '_' savetime '.m'])




end

