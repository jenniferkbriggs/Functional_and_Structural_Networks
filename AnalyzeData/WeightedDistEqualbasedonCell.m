function [Length, syncD_all,nonsyncD_all] = WeightedDistEqualbasedonCell(adj,adjw,tit)
%%This function calculates the shortest weighted distance btw cells that
%%are equal distances away.
% INPUTS:
% - adj: adjacency matrix of the network
% - adjw: weighted structural network  (can be created using weighted net)
% - tit: title for saving

% OUTPUTS: 
% - Length: number of cells separated for corresponding syncD_all and
% nonsyncD_all
% - syncD_all: weighted distances for every synchronzied cell pair
% corresponding to array Length
% - nonsyncD_all: weighted distances for every nonsynchronized cell pair
% corresponding to array Length
numcells = length(adj);
Z = zeros(113, numcells);
ct = 1
AA = sparse(adjw);
GJ = graph(AA, 'upper'); %creates a graph variable using the weighted matrix

for i = 1:numcells %go through each cell and find it's pairs
    if mod(i,200) == 0
        disp(['Cell Number: ' num2str(i)]) %just a way to keep track of progress
    end

    %find connections to cell (i) in the structural and functional
    %networks.
    gj = find(adjw(i,:)); %Weighted structural network
    sync = find(adj(i,:)); %Functional network

    %% Analyze cells longer than 1 GJ away

    %calculate the weighted shortest path between synchronized cell pairs
    for j = 1:(length(sync))
        [Sync(j).data,syncD(j)] = shortestpath(GJ, i, sync(j));
    end

    %calculate the weighted shortest path between nonsynchronized cell
    %pairs
    nonsync = setdiff([1:numcells],sync); %find nonsynchronized pairs
    nonsync(find(nonsync == i))= [];
    for k = 1:length(nonsync)
        [path, d] = shortestpath(GJ, i, nonsync(k));
        Nonsync(k).data = path;
        nonsyncD(k) = d;
    end
    try
        %itterate through all integer distances (e.g. go through 4 (or some number)
        % cells to get from cell (i) to it's synchronized pair). 
        for jj = 1:length(Sync)
            D =[];
            Reflength = length(Sync(jj).data);
            firstref = 1;
            for kk = 1:length(Nonsync)
                if length(Nonsync(kk).data) == Reflength;
                    nonsyncD_all(ct) = nonsyncD(kk); %find all cell pairs RefLength distance from cell (i) that aren't synchronized to cell (i) 
                    if firstref == 1
                        syncD_all(ct) = syncD(jj);
                        firstref = 0; %comment out this line if you want to match non-sync to sync cells by duplicating sync
                    else
                        syncD_all(ct) = NaN;
                    end
                    Length(ct) = Reflength;
                    ct = ct+1;
                end
            end
        end

    end

    clearvars Sync syncD path d Nonsync nonsyncD D
end
try
    save(tit, 'Length', 'syncD_all', 'nonsyncD_all', '-append')
end
end
