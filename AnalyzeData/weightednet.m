function [dist,adjw] = weightednet(gj2adj, gjconduct,example)
%calculates the shortest path of a weighted network
% gj2adj: Adjacency Matrix
% gjconduct: Weighting variable
% example: binary variable which allow you to visualize the analysis based
% on a simple adjacency matrix

%OUTPUTS: 
% -- dist: weighted distances between cell
% -- adjw: adjacency matrix weighted by gjconduct

% Jennifer Briggs 2021
    if example == 1
        gj2adj = [0 1 1 0; 1 0 0 0; 1 0 0 1; 0 0 1 0],
        gjconduct = [1 2 3 2]
    end
    w2 = (gjconduct+gjconduct')/2;
    w = 1./w2;
    adjw = gj2adj.*w;    
    adjw(isnan(adjw))=0;
    A = sparse(adjw);
    dist = distances(graph(A, 'upper'));
    
    if example == 1
       view(biograph(A,[],'ShowWeights','on'))
    end

end
