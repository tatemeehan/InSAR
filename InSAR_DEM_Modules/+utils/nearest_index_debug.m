function [idxMap, dMinMap] = nearest_index_debug(X,Y,Z,P)
    [ny,nx]=size(X); idxMap=nan(ny,nx); dMinMap=nan(ny,nx);
    kdt = KDTreeSearcher(P);
    % sample sparsely to keep fast
    ry = 1:ny; rx = 1:nx; [RR,CC]=ndgrid(ry,rx);
    G = [X(sub2ind([ny,nx],RR(:),CC(:))) Y(sub2ind([ny,nx],RR(:),CC(:))) Z(sub2ind([ny,nx],RR(:),CC(:)))];
    idx = knnsearch(kdt, G);
    D = vecnorm(G - P(idx,:), 2, 2);
    idxMap(sub2ind([ny,nx],RR(:),CC(:))) = idx;
    dMinMap(sub2ind([ny,nx],RR(:),CC(:))) = D;
end