function [L,Lrand,Eglob, Erand, E_loc, nopath, dist,C_alt, C_altrand] = graphProp(adj,x,y,SparseAdjustment)
    numcell = length(adj);
    %Make a random graph


   if ~isempty(x)
   for i = 1:numcell
       for j = 1:numcell
        N1(i,j) = sqrt((x(i)-x(j)).^2+(y(i)-y(j)).^2);
       end
   end
   N11 = reshape(N1,1,numcell^2);
   N = normalize(N11);
   N = reshape(N1, numcell, numcell);
   N = N - diag(diag(N));
   %Neighbor is any cell with an SD < -1
   N = normalize(N1) < -1;

    Iij = N1.*adj;
    dist = mean(nonzeros(Iij),'omitnan')/(mean(nonzeros(N1),'omitnan'));
   else 
       dist = NaN
   end
   
%% Create Random Graph 

  b = zeros(numcell,numcell); %generate an graph with nodes = # cells
  kavg = mean(nansum(adj));
  p = (1000*kavg/numcell);            %p is probability that the cells are connected
 
  for i = 1:numcell
      %i
      for j = i+1:numcell
          %j
          chance = randi(1000);
          if mod(j,2) == 0
              pt = floor(p);
          else
              pt = ceil(p);
          end
          if chance <= pt %if a random number is fits the probability then connect the nodes
              b(i,j) = 1;
              b(j,i) = 1;
          end
      end
  end
            
%% Clustering coefficient.
        % subgraphs of G induced by the neighbors of each vertex
        [MClosed,kClosed,MOpen,kOpen] = subgraphs(adj);
        % local clustering coefficient in each subgraph
        [~,CLocOpen] = deal(zeros(numcell,1));
        [~,CLocClosed] = deal(zeros(numcell,1));

        for i = 1:numcell
            CLocClosed(i) = sum(MClosed{i}(:))/...
                (numel(kClosed{i})*(numel(kClosed{i})-1));
            CLocOpen(i) = sum(MOpen{i}(:))/...
                (numel(kOpen{i})*(numel(kOpen{i})-1));
        end
        % clustering coefficients
%         CClosed = mean(~isnan(CLocClosed));
        CLocOpen(isnan(CLocOpen)) = 0;
        C_alt = mean(CLocOpen);
       % C_alt = mean(CLocOpen, 'omitNaN');
       % C_alt = mean(CLocClosed, 'omitNaN');

  
%% Clustering of Random Graph
        [MClosed_r,kClosed_r,MOpen_r,kOpen_r] = subgraphs(b);
        % local clustering coefficient in each subgraph
        [~,CLocOpen_r] = deal(zeros(numcell,1));
        [~,ClocClosed_r] = deal(zeros(numcell,1));

        for i = 1:numcell
            CLocClosed_r(i) = sum(MClosed_r{i}(:))/...
                (numel(kClosed_r{i})*(numel(kClosed_r{i})-1));
            CLocOpen_r(i) = sum(MOpen_r{i}(:))/...
                (numel(kOpen_r{i})*(numel(kOpen_r{i})-1));
        end
        % clustering coefficients
%         CClosed = mean(~isnan(CLocClosed));
        CLocOpen_r(isnan(CLocOpen_r)) = 0;
        C_altrand = mean(CLocOpen_r);
       % C_altrand = mean(CLocOpen_r,'omitnan');
     %       C_altrand = mean(CLocClosed,'omitnan');
   links = sum(sum(adj))/2;
   linksrand=sum(sum(b))/2;
    for i = 1:numcell
       krand(i) =sum(b(i,:));
    end
    krand = mean(krand);
    
    if abs(kavg-krand)>5;
        disp('Error in random graph matching')
        return
    end
%% Efficiency and Path Lenght    
    A = sparse(adj);
    if SparseAdjustment 
        
    D = graphallshortestpaths(A);
    nopath = (sum(sum(isinf(D))))/(length(D)*length(D-1));%This the the number of pairs of cells with no path
    %Calculate Global Efficiency
    Eglob = (nansum(nansum(1./(D+eye(numcell)))) - numcell)/(numcell*(numcell-1));

    %Calculate Local Efficiency
    [~,ELocOpen] = deal(zeros(numcell,1));

    for i = 1:numcell
       
        if ~isempty(MOpen{i}(:)) & length(MOpen{i}(:))>1
            
            Ab = sparse(MOpen{i});
        LL = graphallshortestpaths(Ab);

        ELocOpen(i) =  (nansum(nansum(1./(LL+eye(length(Ab))))) - length(Ab))...
            /(length(Ab)*(length(Ab)-1));
        
%         sum(Ab)/...
%             (numel(kOpen{i})*(numel(kOpen{i})-1));
        else
         ELocOpen(i) = 0;
        end
    end
        % clustering coefficients
%         CClosed = mean(~isnan(CLocClosed));
    E_loc = mean(ELocOpen);
        
        
    % characteristic path length
    D(isinf(D))= numcell+1;
    L = mean(D(:))/(numcell);
    % global efficiency
% 
    B = sparse(b);
    D = graphallshortestpaths(B);
    %nopath = (sum(sum(isinf(D))))/(length(D)*length(D-1)/2);%This the the number of pairs of cells with no path
    Erand = (nansum(nansum(1./(D+eye(numcell)))) - numcell)/(numcell*(numcell-1));
    D(isinf(D))= numcell+1;
    % global efficiency
    Lrand = mean(D(:))/(numcell);
    else
    
    D = graphallshortestpaths(A);
    nopath = (sum(sum(isinf(D))))/(length(D)*length(D-1));%This the the number of pairs of cells with no path
    Eglob = (nansum(nansum(1./(D+eye(numcell)))) - numcell)/(numcell*(numcell-1));

        for i = 1:numcell
     if ~isempty(MOpen{i}(:)) & length(MOpen{i}(:))>1
            
        Ab = sparse(MOpen{i});
        LL = graphallshortestpaths(Ab);

        ELocOpen(i) =  (nansum(nansum(1./(LL+eye(length(Ab))))) - length(Ab))...
            /(length(Ab)*(length(Ab)-1));
        
%         sum(Ab)/...
%             (numel(kOpen{i})*(numel(kOpen{i})-1));
        else
         ELocOpen(i) = 0;
        end
        end
        % clustering coefficients
%         CClosed = mean(~isnan(CLocClosed));
    E_loc = mean(ELocOpen);
    D(isinf(D))= NaN;
    % characteristic path length
    L = nansum(D(:))/(numcell*(numcell-1));
    % global efficiency
% 
    B = sparse(b);
    D = graphallshortestpaths(B);
    %nopath = (sum(sum(isinf(D))))/(length(D)*length(D-1)/2);%This the the number of pairs of cells with no path
    D(isinf(D))= NaN;
    % global efficiency
    Lrand = nansum(D(:))/(numcell*(numcell-1));
    Erand = (nansum(nansum(1./(D+eye(numcell)))) - numcell)/(numcell*(numcell-1));
    end
    
    
% %
%     figure
%     subplot(2,1,1)
%     imagesc(adj);
%     title('Adjacency Matrix for Islet')
%     subplot(2,1,2)
%     imagesc(b);
%     title('Random Graph with same nodes, links, and avg k')

end
function [MClosed,kClosed,MOpen,kOpen] = subgraphs(A)
% subgraphs: compute adjacency matrices for each vertex in a graph
% usage: [MClosed,kClosed,MOpen,kOpen] = subgraphs(A);
%
% arguments: 
%   A - (nxn) adjacency matrix of a graph G
%
%   MClosed, MOpen - (nx1) cell arrays containing adjacency matrices of the 
%       subgraphs corresponding to neighbors of each vertex. For example, 
%       MClosed{j} is the adjacency matrix of the subgraph of G 
%       corresponding to the closed neighborhood of the jth vertex of G, 
%       and kClosed{j} is the list of vertices of G that are in the 
%       subgraph (and represent the corresponding rows/columns of 
%       MClosed{j})
%       
% author: Nathan D. Cahill
% email: nathan.cahill@rit.edu
% date: 10 Apr 2014
% number of vertices
n = size(A,1);
% initialize indices of neighboring vertices, and adjacency matrices
[kClosed,kOpen] = deal(cell(n,1));
[MClosed,MOpen] = deal(cell(n,1));
% loop through each vertex, determining its neighbors
for i=1:n
    
    % find indices of neighbors of vertex v_i
    k1 = find(A(i,:)); % vertices with edge beginning at v_i
    k2 = find(A(:,i)); % vertices with edge ending at v_i
    kOpen{i} = union(k1,k2); % vertices in neighborhood of v_i
     kClosed{i} = union(kOpen{i},i);
    
    % extract submatrix describing adjacency of neighbors of v_i
    MOpen{i} = A(kOpen{i},kOpen{i});
     MClosed{i} = A(kClosed{i},kClosed{i});
    
end
end