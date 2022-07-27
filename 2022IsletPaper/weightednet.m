function [dist,adjw] = weightednet(gj2adj, gjconduct,gjconn)
%calculates the shortest path of a weighted network
% gj2adj: Adjacency Matrix
% gjconduct: Weighting variable
% gjconn: binary variable which allow you to visualize the analysis based
% on a simple adjacency matrix
% Jennifer Briggs 2021
    if gjconn == 1
        gj2adj = [0 1 1 0; 1 0 0 0; 1 0 0 1; 0 0 1 0],
        gjconduct = [1 2 3 2]
    end
    w2 = (gjconduct+gjconduct')/2;
    w = 1./w2;
    adjw = gj2adj.*w;        
    A = sparse(adjw);
    dist = graphallshortestpaths(A, 'Directed', 'False');
    
    if gjconn == 1
       view(biograph(A,[],'ShowWeights','on'))
    end

end
