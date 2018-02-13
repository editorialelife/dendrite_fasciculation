function [bwPath3D,pt1,pt2] = brightestPathBetweenTwoPts3D(im, xyPxSize, zPxSize, varargin)
%BRIGHTESTPATHBETWEENTWOPTS extracts highest intensity path between pts.

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('im');
ip.addRequired('xyPxSize');
ip.addRequired('zPxSize');
ip.addOptional('pt1',0);
ip.addOptional('pt2',0);
ip.addOptional('filterNoise', 1);
ip.addOptional('filterBackground', 1);
ip.addOptional('shrinkFactor', 1);
% ip.addOptional('mask', -1);
ip.parse(im, xyPxSize, zPxSize, varargin{:});
pt1 = ip.Results.pt1;
pt2 = ip.Results.pt2;
filterNoise = ip.Results.filterNoise;
filterBackground = ip.Results.filterBackground;
shrinkFactor = ip.Results.shrinkFactor;
% mask = ip.Results.mask;

% check the dimension of the input dataset
nd = numel(size(im));
if nd ~= 3
    error(message('This function only takes in 3D image data'));
end   

% make sure coordinates of both seeding pts are within the limit of the
% dataset
[ySize, xSize, zSize] = size(im);
if numel(pt1) == 3
    if pt1(1) < 1 || pt1(1) > ySize || pt1(2) < 1 || pt1(2) > xSize...
            || pt2(1) < 1 || pt2(1) > ySize || pt2(2) < 1 || pt2(2) > xSize ...
            || pt1(3) < 1 || pt1(3) > zSize || pt1(3) < 1 || pt1(3) > zSize
        error(message('Seeding point coordinate out of bounds'));
    end
end

% convert image to double
im = im2double(im);

% filter noise for each slice
if filterNoise
    imF = zeros(size(im));
    for z = 1:zSize
        imF(:,:,z) = filterGauss2D(im(:,:,z),1);
    end
else
    imF = im;
end

% filter background for each slice
imB = zeros(size(im));
if filterBackground
    for z = 1:zSize
        imB(:,:,z) = filterGauss2D(im(:,:,z),10);
    end
end
imCorrected = imF - imB;

%% first detect the brightest path in z-projection
% z-projection
imPrj = max(imCorrected,[],3);
% perform correction on the projection image again
if filterNoise
    imPrj = filterGauss2D(imPrj,1);
end
if filterBackground
    imPrj = imPrj - filterGauss2D(imPrj,10);
end
% normalize the intensity range to [0,1]
imPrj = (imPrj - min(imPrj(:)))/(max(imPrj(:)) - min(imPrj(:)));

if shrinkFactor < 1
    imPrjS=imresize(imPrj,shrinkFactor);
    [ySizeS, xSizeS, ~] = size(imPrjS);
    pt1S(1) = min(max(round(pt1(1) * shrinkFactor),1), ySizeS);
    pt1S(2) = min(max(round(pt1(2) * shrinkFactor),1), xSizeS);
    pt2S(1) = min(max(round(pt2(1) * shrinkFactor),1), ySizeS);
    pt2S(2) = min(max(round(pt2(2) * shrinkFactor),1), xSizeS);
    pathXYS = findBrightestPathGraph2D(imPrjS, pt1S, pt2S, true(size(imPrjS)));
    tmp=false(size(imPrjS));
    tmp(pathXYS)=true;
    mask=imdilate(imresize(tmp,size(imPrj)),strel('disk',ceil(1/shrinkFactor)+1));
end


% % threshold z-projection image to create a 2D mask
% if mask == -1
%     [~,th] = cutFirstHistMode(imPrj,0);
%     mask = imPrj > th;
%     % refine the threshold so the two points are enclosed in one mask
%     l = bwlabel(mask);
%     idx1 = l(pt1(1),pt1(2));
%     idx2 = l(pt2(1),pt2(2));
%     if idx1 == 0 || idx2 == 0 || idx1 ~= idx2
%         thresholdFactor = 0.9;
%         while idx1 == 0 || idx2 == 0 || idx1 ~= idx2
%             th = th * thresholdFactor;
%             mask = imPrj > th;
%             l = bwlabel(mask);
%             idx1 = l(pt1(1),pt1(2));
%             idx2 = l(pt2(1),pt2(2));
%         end
%     else
%         thresholdFactor = 1.1;
%         while idx1 > 0 && idx1 == idx2
%             th = th * thresholdFactor;
%             mask = imPrj > th;
%             l = bwlabel(mask);
%             idx1 = l(pt1(1),pt1(2));
%             idx2 = l(pt2(1),pt2(2));
%         end
%         th = th / thresholdFactor;
%     end
%     % only keep the patch of mask that enclose both tips
%     mask = imPrj > th;
%     l = bwlabel(mask);
%     idx = l(pt1(1),pt1(2));
%     mask = l == idx;
%     mask = bwmorph(mask,'diag');
% end
% % update figure
% centDisp = false(size(imPrj));
% centDisp(pt1(1),pt1(2)) = true;
% centDisp(pt2(1),pt2(2)) = true;
% centDisp=imdilate(centDisp,strel('disk',2));
% imDisp = cat(3,max(imPrj,0.5*double(mask))-centDisp,max(imPrj,centDisp),max(imPrj,centDisp));
% imDisp(imDisp < 0) = 0;
% figure,imshow(imDisp);

% call function to find the brightest path in the x-y plane from the
% projection image
pathXY = findBrightestPathGraph2D(imPrj, pt1, pt2, mask);
bwPathXY = false(size(imPrj));
bwPathXY(pathXY) = true;
bwPathXY = bwmorph(bwPathXY,'skel','inf');
% extract sequential pixels based on geodesic distance
bwPt1XY = false(size(bwPathXY));
bwPt1XY(pt1(1),pt1(2))=true;
gdDist = bwdistgeodesic(bwPathXY,bwPt1XY);
ptsInd = find(bwPathXY);
distVals = gdDist(ptsInd);
tmp=cat(2,distVals,ptsInd);
tmp=sortrows(tmp,1);
pathXY = tmp(:,2);
nPx = sum(bwPathXY(:));

%% find the brightest path in Z
% extract the trace from each z-slice and combine into a 2D image
zCut = zeros(zSize,nPx);
for z = 1:zSize
    tmp = imCorrected(:,:,z);
    zCut(z,:) = tmp(pathXY);
end
zCut = zCut-min(zCut(:));
% interpolate in z-axis
zMax = zPxSize*zSize;
zCut=interp1(zPxSize:zPxSize:zMax,zCut,zPxSize:xyPxSize:zMax,'linear');
% figure out the start and end z position
if numel(pt1) == 2
    [~,pt1(3)] = max(zCut(:,1));
    [~,pt2(3)] = max(zCut(:,nPx));
end
pt1Z = [pt1(3),1];
pt2Z = [pt2(3),nPx];
% call function to find the brightest path in Z
pathZ = findBrightestPathGraph2D(zCut,pt1Z,pt2Z);
bwPathZ = false(size(zCut));
bwPathZ(pathZ) = true;
bwPathZ = bwmorph(bwPathZ,'skel','inf');
% extract sequential pixels based on geodesic distance
bwPt1Z = false(size(bwPathZ));
bwPt1Z(pt1Z(1),pt1Z(2))=true;
gdDist = bwdistgeodesic(bwPathZ,bwPt1Z);
ptsInd = find(bwPathZ);
distVals = gdDist(ptsInd);
tmp=cat(2,distVals,ptsInd);
tmp=sortrows(tmp,1);
pathZ = tmp(:,2);

% pathZ = find(bwPathZ);

%% reconstruct the brightest path in 3D
[zPos, xyInd] = ind2sub(size(zCut), pathZ);
pathXY2 = zeros(numel(xyInd),1);
for i = 1 : numel(xyInd)
    pathXY2(i) = pathXY(xyInd(i));
end
[yPos, xPos] = ind2sub(size(imPrj), pathXY2);
bwPath3D = false(size(im,1),size(im,2),size(zCut,1));
for i = 1:numel(xyInd)
    bwPath3D(yPos(i),xPos(i),zPos(i)) = true;
end

% path3D = sub2ind([size(im,1),size(im,2),size(zCut,1)], yPos, xPos, zPos);
% bwPath3D = false(size(im,1),size(im,2),size(zCut,1));
% bwPath3D(path3D) = true;

end
