function path = findBrightestPathGraph2D(im, pt1, pt2, mask)

% check input
if nargin < 4
    mask = true(size(im));
end

% check the dimension of the input dataset
nd = numel(size(im));
if nd ~= 2
    error(message('This function only takes in 2D image data'));
end

% loop through all the pixels to find adjcent pairs, the weight of the edge
% is the average intensity of the two pixels
[ySize, xSize] = size(im);
nNode = sum(mask(:));
pxIndList = find(mask);
% nNode = ySize * xSize * zSize;
edge = [];
w = [];
n = 0;
for i = 1 : nNode-1
    [y, x] = ind2sub([ySize, xSize], pxIndList(i));
    % x + 1
    if x < xSize
        newPxInd = sub2ind([ySize, xSize], y, x+1);
        idx = find(pxIndList == newPxInd);
        if ~isempty(idx)
            n = n + 1;
            edge(n,1) = i;
            edge(n,2) = idx;
            w(n) = (im(y,x) + im(y,x+1))/2;
        end
    end
    % y + 1
    if y < ySize
    newPxInd = sub2ind([ySize, xSize], y+1, x);
    idx = find(pxIndList == newPxInd);
    if ~isempty(idx)
        n = n + 1;
        edge(n,1) = i;
        edge(n,2) = idx;
        w(n) = (im(y,x) + im(y+1,x))/2;
    end
    end
end

% to use the shortest path function, the weight is the 1/pixel intensity in
% order to find the brightest pixels along the path
w = (1 ./ w)';
% w = - w';

% contrust the adjacency matrix
adjMat = sparse(vertcat(edge(:,1), edge(:,2)), vertcat(edge(:,2), edge(:,1)), repmat(w, [2 1]));

% convert the coordinate of the two seeded pts into indices
pt1Ind = sub2ind([ySize, xSize], pt1(1), pt1(2));
pt2Ind = sub2ind([ySize, xSize], pt2(1), pt2(2));
pt1IndOrder = find(pxIndList == pt1Ind);
pt2IndOrder = find(pxIndList == pt2Ind);

% call shortest path function to find the brightest path
[~, path, ~] = graphshortestpath(adjMat, pt1IndOrder, pt2IndOrder, 'Directed', false);
nPts = numel(path);
for n = 1:nPts
    path(n) = pxIndList(path(n));
end

end