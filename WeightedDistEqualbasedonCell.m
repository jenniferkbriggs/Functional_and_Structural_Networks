function [Length, syncD_all,nonsyncD_all] = WeightedDistEqualbasedonCell(adj,adjw,tit)
%%This function calculates the shortest weighted distance btw cells with
%%equal distances away.
% figure
% xlabel('Cells')
% ylabel('Average Z')
numcells = length(adj);
Z = zeros(113, numcells);
ct = 1
AA = sparse(adjw);
GJ = graph(AA, 'upper');

for i = 1:numcells
    if mod(i,200) == 0
        disp(['Cell Number: ' num2str(i)])
    end
    gj = find(adjw(i,:)); %Weighted structural network
    sync = find(adj(i,:)); %Functional network

    %% Analyze cells longer than 1 GJ away
    for j = 1:(length(sync))
        [Sync(j).data,syncD(j)] = shortestpath(GJ, i, sync(j));
    end

    nonsync = setdiff([1:numcells],sync);
    nonsync(find(nonsync == i))= [];
    for k = 1:length(nonsync)
    [path, d] = shortestpath(GJ, i, nonsync(k));
    Nonsync(k).data = path;
    nonsyncD(k) = d;
    end
    try
        for jj = 1:length(Sync)
            D =[];
            Reflength = length(Sync(jj).data);
            firstref = 1;
            for kk = 1:length(Nonsync)
                if length(Nonsync(kk).data) == Reflength;
                    nonsyncD_all(ct) = nonsyncD(kk);
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
        %if i == 1
        %    save('ShortPathEqualDist2', 'Length', 'syncD_all', 'nonsyncD_all')
        %    elseif mod(i,100) == 0
        %        save('ShortPathEqualDist2', 'Length', 'syncD_all', 'nonsyncD_all', '-append')
        %end
    end
    % hold on
    % plot(i,mean(nonzeros(Z),'omitnan'), 'o')
    % drawnow
    clearvars Sync syncD path d Nonsync nonsyncD D
    end
    save(tit, 'Length', 'syncD_all', 'nonsyncD_all', '-append')
end
