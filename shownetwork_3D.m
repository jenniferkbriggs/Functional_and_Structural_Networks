function out = shownetwork_3D(adj, pos, both, gj2adj, tit, savename,hubby)
%function to view 3D network
%Jennifer Briggs 2020
%Inputs:
%adj = adjacency matrix
%pos = 3D position matrix
%both = graph types
savetime =  datestr(datetime('today'),'yyyymmdd');
if 0

importdata([pwd '\' savename '_' savetime '.m'])
names = fieldnames(ans)
indx = strfind(names, 'ans')
for i = 1:length(names)
    if indx{i} ~= 1
    eval([names{i} '=ans.' names{i}]);
    end
end
end

prompt = 'Plot Hub [H], Random [R], High GJ [GJ], Low connections [L]'
str = input(prompt, 's');
AdjacencyGraph = graph(adj, 'upper');
x = pos(:,1);
y = pos(:,2);
z = pos(:,3);
numcell = length(adj)
if str == 'H'
repcell = randi([1 length(hubby)],1,1);
repcell = hubby(repcell)
elseif str == 'GJ'
[~, repcell] = max(sum(gj2adj))
repcell = find(sum(gj2adj) == round(mean(sum(gj2adj)))+3)
repcell = repcell(4) ;
elseif str == 'R'
repcell = round(rand().*numcell);
elseif str == 'L'
s = sum(adj, 'omitnan');
s(s == 0) = NaN;
[~, repcell] = min(s)
end


cellmat = zeros(numcell,numcell);
cellmat(repcell,:) = 1;
cellmat(:,repcell) = 1;
both2 = both.*cellmat;
gj2adj2 = gj2adj.*cellmat;
adj2 = adj.*cellmat;

g = graph(adj2+gj2adj2+both2,'upper');
g1 = graph(adj2,'upper');

g1c = repmat([0 0 0],numcell,1); %RED
g1c((adj2 (:,repcell) == 1),1) = 0.01%[0.6350];
g1c((adj2 (:,repcell) == 1),2) = 0.01%0.0780;
g1c((adj2 (:,repcell) == 1),3)= 1%0.1840;

g3c = repmat([[0 0 0]],numcell,1); %BLACK
g3c((both2 (:,repcell) == 1),1) = [0];
g3c((both2 (:,repcell) == 1),2) = 0;
g3c((both2 (:,repcell) == 1),3)= 0;

g2c = repmat([[0 0 0]],numcell,1); %GREEN
g2c((gj2adj2 (:,repcell) == 1),1) =  1%0.3650;
g2c((gj2adj2 (:,repcell) == 1),2) = .01%0.9220;
g2c((gj2adj2 (:,repcell) == 1),3)= .01%0.8160;

indx2 = find(gj2adj2 (:,repcell) == 1)
indx1 = find(adj2 (:,repcell) == 1)
indxboth = intersect(indx2,indx1)
% for lll 1:length(indx2)
    
gc = g1c+g2c;
gc(gc == 0) = .9;
gc(gc == 1) = 0;

g1c(g1c == 0) = .9;
g2c(g2c == 0) = .9;
g3c(g3c == 0) = .9;
g3c(g3c == 1) = 0;

% cell1 = find(g1c(:,1) ~= 0.6350);
% cell2 = find(g2c(:,1) ~= 0.3650);
cell1 = find(g1c(:,1) ~= 0.01);
cell2 = find(g2c(:,1) ~= 1);

cell2 = intersect(cell1,cell2);
cell2 = sort([cell2; indxboth])

g2 = graph(gj2adj2(cell1,cell1),'upper');
g3 = graph(both2(cell2,cell2),'upper');

Edgec = repmat([.175 .54 .60],length(table2array(g.Edges)),1);
e = table2array(g.Edges);
for i=1:length(Edgec)
    if find(indxboth == [e(i,1)  e(i,2)]);
        Edgec(i,:) = [ 0 0 0];
    elseif find(indx1 == [e(i,1)  e(i,2)]);
        Edgec(i,:) = [0.01 0.01 1]; %[0.6350 0.0780 0.1840]
    elseif find(indx2 == [e(i,1)  e(i,2)]);
        Edgec(i,:) = [1 0.01 0.01];%[0.3650 0.9220 0.8160]
    end
end
% 
% Edgec1 = repmat([0.6350 0.0780 0.1840],size(table2array(g1.Edges),1),1);
% Edgec2 = repmat([.175 .54 .60],size(table2array(g2.Edges),1),1);

Edgec1 = repmat([0.01 0.01 1],size(table2array(g1.Edges),1),1);
Edgec2 = repmat([1 0.01 0.01],size(table2array(g2.Edges),1),1);
Edgec3 = repmat([0 0 0],size(table2array(g3.Edges),1),1);


g3c(repcell,:) = [0 0 0];
for ll = 1:length(indxboth)
g3c(indxboth(ll), :) = [0 0 0];
end
Msize = ones(numcell, 1).*5;
Msize(gc(:,1) ~= .9) = 10;
Msize(repcell) = 15;
% Nnames = strings(numcell,1)
% % Nnames(gc(:,1) == 0) = 'Both'
% % Nnames(gc(:,1) == 0.3650) = 'GJ'
% % Nnames(gc(:,1) == 0.6350) = 'Sync'
% Nnames(indxboth(1)) = 'Both';
% Nnames(indx1(3)) = 'Sync';
% Nnames(indx2(5)) = 'GJ';
% 


figure, plot(g1,'Xdata',x,'YData',y,'ZData', z,'EdgeColor', Edgec1, ...
    'NodeColor',g1c, 'LineWidth', 2, 'MarkerSize', Msize)
hold on, plot(g2,'Xdata',x(cell1),'YData',y(cell1),'ZData', z(cell1),'EdgeColor', Edgec2, ...
    'NodeColor',g2c(cell1,:), 'LineWidth', 2, 'MarkerSize', Msize(cell1))
hold on, plot(g3,'Xdata',x(cell2),'YData',y(cell2),'ZData', z(cell2),'EdgeColor', Edgec3, ...
    'NodeColor',g3c(cell2,:), 'LineWidth', 2, 'MarkerSize', Msize(cell2))
% 
% figure, plot(g, 'Xdata',x,'YData',y,'ZData', z,'EdgeColor', Edgec, ...
%     'NodeColor',g1c, 'LineWidth', 2, 'MarkerSize', Msize, 'NodeLabel', Nnames)
% hold on
% plot(g2, 'Xdata',x,'YData',y,'ZData', z)
set(gca,'ytick',[])
set(gca,'xtick',[])
set(gca,'ztick',[])
legend( '', '', '','Synchronized', 'Gap Junction', 'Both', 'FontSize', 20)
title(tit,'FontSize', 30)
saveas(gcf, [pwd '\' savename '_' savetime '.fig'])
saveas(gcf, [pwd '\' savename '_' savetime '.png'])
save([pwd '\' savename '_' savetime '.m'])
end