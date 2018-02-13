function neuralProcessesSpacialRelationshipBatch_TIF(dataPath,varargin)

%% input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPath');
ip.addParamValue('xyPxSize',  0.11, @(x) isnumeric(x));
ip.addParamValue('zPxSize', 0.35, @(x) isnumeric(x));
ip.parse(dataPath, varargin{:});
xyPxSize = ip.Results.xyPxSize;
zPxSize = ip.Results.zPxSize;

if ~exist(dataPath,'dir')
    dataPath = uigetdir(pwd,'Please select folder that contains input images in the inputImages subfolder.');
end
dataPath = [dataPath filesep];

% input images path
inputPath = [dataPath 'inputImages' filesep];
if ~exist(inputPath,'dir')
    error('No inputImages subfolder found!');
end

% output path
outputPath = [dataPath 'results' filesep];
if ~exist(outputPath,'dir');
    mkdir(outputPath);
end

% load in input images list
fileList1 = dir([inputPath '*.tif']);
fileList2 = dir([inputPath '*.TIF']);
fileList = cat(1,fileList1,fileList2);
nDatasets = numel(fileList);
display(['Total of ' num2str(nDatasets) ' to be processed'])

if nDatasets == 0
    warning('No images to be processed!');
else
    %% looping over all the dataset to do the initial seeding of tip points
    % create folder to store seed points
    tempPath = [dataPath 'seedPts' filesep];
    if ~exist(tempPath,'dir')
        mkdir(tempPath);
    end
    % loop over all datasets to do initial seeding
    tempFilename = cell(nDatasets);
    for n = 1:nDatasets
        % load in three channels
        filename = [inputPath fileList(n).name];
        [blue,green,red] = loadImagesNeuralProcessesSpacialRelationship_TIF(filename);
        [ySize, xSize, zSize] = size(blue);
        
        % no cropping in batch mode
        bCropBound = [1,zSize];
        gCropBound = [1,zSize];
        rCropBound = [1,zSize];
        cropBound(1) = min(min(bCropBound(1),gCropBound(1)),rCropBound(1));
        cropBound(2) = max(max(bCropBound(2),gCropBound(2)),rCropBound(2));
        blue = blue(:,:,bCropBound(1):bCropBound(2));
        green = green(:,:,gCropBound(1):gCropBound(2));
        red = red(:,:,rCropBound(1):rCropBound(2));
        
        % tip and base points for each channel
        [bPt1,bPt2] = pickPointsFromProjectionImage(blue);
        [gPt1,gPt2] = pickPointsFromProjectionImage(green);
        [rPt1,rPt2] = pickPointsFromProjectionImage(red);
        
        % save everything to file
        tmp = fileList(n).name;
        tempFilename{n} = [tmp(1:end-3) '_seed_points.mat'];
        save([tempPath tempFilename{n}],'bPt1','bPt2','gPt1','gPt2','rPt1','rPt2');
    end
    
    %% proceed with all detection steps
    display('Starting batch processing data');
    for n = 1:nDatasets
        display(['Processing dataset: ' fileList(n).name]);
        %% load in data
        % load in images
        filename = [inputPath fileList(n).name];
        [blue,green,red] = loadImagesNeuralProcessesSpacialRelationship_TIF(filename);
        [ySize, xSize, zSize] = size(blue);
        
        % no cropping in batch mode
        bCropBound = [1,zSize];
        gCropBound = [1,zSize];
        rCropBound = [1,zSize];
        cropBound(1) = min(min(bCropBound(1),gCropBound(1)),rCropBound(1));
        cropBound(2) = max(max(bCropBound(2),gCropBound(2)),rCropBound(2));
        blueOrig = blue(:,:,cropBound(1):cropBound(2));
        greenOrig = green(:,:,cropBound(1):cropBound(2));
        redOrig = red(:,:,cropBound(1):cropBound(2));
        blue = blue(:,:,bCropBound(1):bCropBound(2));
        green = green(:,:,gCropBound(1):gCropBound(2));
        red = red(:,:,rCropBound(1):rCropBound(2));
        
        % load in data for the current set
        load([tempPath tempFilename{n}]);
        
        %% extract traces from 3 channels
        % call function to extract the bright trace from 3 channels
        tic; [bwBPath3D,bPt1,bPt2] = brightestPathBetweenTwoPts3D(blue,xyPxSize,zPxSize,bPt1,bPt2,0,0,0.25); t = toc;
        display(['Finished extracting trace from the blue channel in ' num2str(t) ' secs'])
        tic; [bwGPath3D,gPt1,gPt2] = brightestPathBetweenTwoPts3D(green,xyPxSize,zPxSize,gPt1,gPt2,0,0,0.25); t = toc;
        display(['Finished extracting trace from the green channel in ' num2str(t) ' secs'])
        tic; [bwRPath3D,rPt1,rPt2] = brightestPathBetweenTwoPts3D(red,xyPxSize,zPxSize,rPt1,rPt2,0,0,0.25); t = toc;
        display(['Finished extracting trace from the red channel in ' num2str(t) ' secs'])
        
        % pad the 3D traces from 3 channels to the same size in Z
        % cropBound(1) = min(min(bCropBound(1),gCropBound(1)),rCropBound(1));
        % cropBound(2) = max(max(bCropBound(2),gCropBound(2)),rCropBound(2));
        tmp = size(blueOrig,3);
        zMax = zPxSize*tmp;
        zSize2 = numel(zPxSize:xyPxSize:zMax);
        tmp = (bCropBound(1)-cropBound(1)+1)*zPxSize;
        tmp = numel(tmp:xyPxSize:zMax);
        bZ1 = zSize2-tmp+1;
        tmp = (gCropBound(1)-cropBound(1)+1)*zPxSize;
        tmp = numel(tmp:xyPxSize:zMax);
        gZ1 = zSize2-tmp+1;
        tmp = (rCropBound(1)-cropBound(1)+1)*zPxSize;
        tmp = numel(tmp:xyPxSize:zMax);
        rZ1 = zSize2-tmp+1;
        % bZ1 = round((bCropBound(1)-cropBound(1))*zPxSize/xyPxSize)+1;
        % gZ1 = round((gCropBound(1)-cropBound(1))*zPxSize/xyPxSize)+1;
        % rZ1 = round((rCropBound(1)-cropBound(1))*zPxSize/xyPxSize)+1;
        tmp = false(ySize,xSize,zSize2);
        tmp(:,:,bZ1:bZ1+size(bwBPath3D,3)-1) = bwBPath3D;
        bwBPath3D = tmp;
        bPt1(3) = bPt1(3) + bZ1 - 1;
        bPt2(3) = bPt2(3) + bZ1 - 1;
        tmp = false(ySize,xSize,zSize2);
        tmp(:,:,gZ1:gZ1+size(bwGPath3D,3)-1) = bwGPath3D;
        bwGPath3D = tmp;
        gPt1(3) = gPt1(3) + gZ1 - 1;
        gPt2(3) = gPt2(3) + gZ1 - 1;
        tmp = false(ySize,xSize,zSize2);
        tmp(:,:,rZ1:rZ1+size(bwRPath3D,3)-1) = bwRPath3D;
        bwRPath3D = tmp;
        rPt1(3) = rPt1(3) + rZ1 - 1;
        rPt2(3) = rPt2(3) + rZ1 - 1;
        
        %% extract the average trace that runs along the center of the 3 traces
        % automatically pick the two tips of the average trace
        bwPaths3D = {bwBPath3D,bwGPath3D,bwRPath3D};
        pt1s = {bPt1,gPt1,rPt1};
        pt2s = {bPt2,gPt2,rPt2};
        [pt1,pt2] = findTipPointsOfAvgTrace(bwPaths3D,pt1s,pt2s);
        
        % call function to extrace the average trace based on summed distances to
        % all 3 traces
        sumDist = bwdist(bwBPath3D);
        sumDist = sumDist+bwdist(bwGPath3D);
        sumDist = sumDist+bwdist(bwRPath3D);
        sumDist = 1./(sumDist+1);
        % bwPath3D = brightestPathBetweenTwoPts3D(sumDist,1,1,pt1,pt2,0,0,combMask);
        tic; bwPath3D = brightestPathBetweenTwoPts3D(sumDist,1,1,pt1,pt2,0,0,0.25); t = toc;
        display(['Finished extracting the average trace in ' num2str(t) ' secs'])
        
        %% walk along the average trace to find calculate the relative relationship of the 3 traces
        % extract the 1st pt
        pt1Z = find(bwPath3D(pt1(1),pt1(2),:));
        pt1Z = pt1Z(1);
        bwPt1 = false(size(bwPath3D));
        bwPt1(pt1(1),pt1(2),pt1Z)=true;
        % compute the geodesic distance of all the pxs to the 1st point and order
        % them
        gdDist = bwdistgeodesic(bwPath3D,bwPt1);
        ptsInd = find(bwPath3D);
        nPts = numel(ptsInd);
        distVals = gdDist(ptsInd);
        [~,orderInd] = sort(distVals,'ascend');
        
        % loop through all the pxs along the average trace
        directDist = zeros(nPts,1);
        geodesicDist = zeros(nPts,1);
        pairDist = zeros(nPts,3);
        midProcInd = zeros(nPts,1);
        colorPath2D = zeros(ySize,xSize);
        for i = 1:nPts
            [y,x,z] = ind2sub(size(bwPath3D),ptsInd(orderInd(i)));
            directDist(i) = sqrt((pt1(1)-y)^2+(pt1(2)-x)^2+(pt1Z-z)^2);
            geodesicDist(i) = gdDist(ptsInd(orderInd(i)));
            bwPtI3D = false(size(bwPath3D));
            bwPtI3D(ptsInd(orderInd(i))) = true;
            
            flag = 0;
            searchR = 10;
            while flag == 0
                ymin = max(y-searchR,1);
                ymax = min(y+searchR,ySize);
                xmin = max(x-searchR,1);
                xmax = min(x+searchR,xSize);
                zmin = max(z-searchR,1);
                zmax = min(z+searchR,zSize2);
                cropB = bwBPath3D(ymin:ymax,xmin:xmax,zmin:zmax);
                cropG = bwGPath3D(ymin:ymax,xmin:xmax,zmin:zmax);
                cropR = bwRPath3D(ymin:ymax,xmin:xmax,zmin:zmax);
                if max(cropB(:)) > 0 && max(cropG(:)) > 0 && max(cropR(:)) > 0
                    crop = bwPtI3D(ymin:ymax,xmin:xmax,zmin:zmax);
                    flag = 1;
                else
                    searchR = searchR+10;
                end
            end
            
            tmp = bwdist(crop);
            tmp2 = tmp.*double(cropB);
            tmp2(~cropB)=max(tmp2(:));
            [~,bInd] = min(tmp2(:));
            [bY,bX,bZ]=ind2sub(size(crop),bInd);
            tmp2 = tmp.*double(cropG);
            tmp2(~cropG)=max(tmp2(:));
            [~,gInd] = min(tmp2(:));
            [gY,gX,gZ]=ind2sub(size(crop),gInd);
            tmp2 = tmp.*double(cropR);
            tmp2(~cropR)=max(tmp2(:));
            [~,rInd] = min(tmp2(:));
            [rY,rX,rZ]=ind2sub(size(crop),rInd);
            rgDist = sqrt((rY-gY)^2+(rX-gX)^2+(rZ-gZ)^2);
            brDist = sqrt((bY-rY)^2+(bX-rX)^2+(bZ-rZ)^2);
            gbDist = sqrt((gY-bY)^2+(gX-bX)^2+(gZ-bZ)^2);
            pairDist(i,:) = [rgDist,brDist,gbDist];
            maxDist = max(max(rgDist,brDist),gbDist);
            if maxDist == rgDist
                midProcInd(i) = 1;  % blue in the middle
            elseif maxDist == brDist
                midProcInd(i) = 2;  % green in the middle
            else
                midProcInd(i) = 3;  % red in the middle
            end
            colorPath2D(y,x) =  midProcInd(i);
        end
        
        %% output
        % make individual output folders
        name = fileList(n).name;
        name = name(1:end-3);
        outPath = [outputPath name '_results' filesep];
        if ~exist(outPath,'dir')
            mkdir(outPath)
        end
        
        % save numerical results to csv file
        outputName = [outPath name '_output.csv'];
        fid1=fopen(outputName,'wt');
        fprintf(fid1,'Pixel_index,Geodesic_dist_end,Direct_dist_end,RG_dist,BR_dist,GB_dist,Middle_process(B=1,G=2,R=3)\n');
        for i = 1:nPts
            fprintf(fid1,'%s,%s,%s,%s,%s,%s,%s\n',num2str(i),num2str(geodesicDist(i)),num2str(directDist(i)),...
                num2str(pairDist(i,1)),num2str(pairDist(i,2)),num2str(pairDist(i,3)),num2str(midProcInd(i)));
        end
        fclose(fid1);
        
        % output images
        colorPath2D = cat(3, colorPath2D == 3, colorPath2D == 2, double(colorPath2D == 1));
        imwrite(colorPath2D, [outPath 'colorPath2D.tif'], 'tif', 'compression', 'none');
        stackWrite(bwPath3D, [outPath 'bwPath3D.tif']);
        stackWrite(bwBPath3D, [outPath 'bwBPath3D.tif']);
        stackWrite(bwGPath3D, [outPath 'bwGPath3D.tif']);
        stackWrite(bwRPath3D, [outPath 'bwRPath3D.tif']);
        
        % X-Y projection overlay of traces
        % blue
        im = imadjust(im2double(max(blue,[],3)));
        bw = max(bwBPath3D,[],3);
        xyPrjOverlay = cat(3,max(im,bw),im-bw,im-bw);
        imwrite(xyPrjOverlay, [outPath 'bPathOverlayXY.tif'], 'tif', 'compression', 'none');
        % green
        im = imadjust(im2double(max(green,[],3)));
        bw = max(bwGPath3D,[],3);
        xyPrjOverlay = cat(3,max(im,bw),im-bw,im-bw);
        imwrite(xyPrjOverlay, [outPath 'gPathOverlayXY.tif'], 'tif', 'compression', 'none');
        % red
        im = imadjust(im2double(max(red,[],3)));
        bw = max(bwRPath3D,[],3);
        xyPrjOverlay = cat(3,max(im,bw),im-bw,im-bw);
        imwrite(xyPrjOverlay, [outPath 'rPathOverlayXY.tif'], 'tif', 'compression', 'none');
        
        % interpolation of original stacks and Y-Z projection overlay of traces
        % blue
        interpStk = permute(interp1(zPxSize:zPxSize:zMax,permute(im2double(blueOrig),[3 1 2]),zPxSize:xyPxSize:zMax,'linear'),[2 3 1]);
        stackWrite(im2uint16(interpStk),[outPath 'blue_interp.tif']);
        im = imadjust(im2double(squeeze(max(interpStk,[],2))));
        bw = squeeze(max(bwBPath3D,[],2));
        yzPrjOverlay = cat(3,max(im,bw),im-bw,im-bw);
        imwrite(yzPrjOverlay, [outPath 'bPathOverlayYZ.tif'], 'tif', 'compression', 'none');
        % green
        interpStk = permute(interp1(zPxSize:zPxSize:zMax,permute(im2double(greenOrig),[3 1 2]),zPxSize:xyPxSize:zMax,'linear'),[2 3 1]);
        stackWrite(im2uint16(interpStk),[outPath 'green_interp.tif']);
        im = imadjust(im2double(squeeze(max(interpStk,[],2))));
        bw = squeeze(max(bwGPath3D,[],2));
        yzPrjOverlay = cat(3,max(im,bw),im-bw,im-bw);
        imwrite(yzPrjOverlay, [outPath 'gPathOverlayYZ.tif'], 'tif', 'compression', 'none');
        % red
        interpStk = permute(interp1(zPxSize:zPxSize:zMax,permute(im2double(redOrig),[3 1 2]),zPxSize:xyPxSize:zMax,'linear'),[2 3 1]);
        stackWrite(im2uint16(interpStk),[outPath 'red_interp.tif']);
        im = imadjust(im2double(squeeze(max(interpStk,[],2))));
        bw = squeeze(max(bwRPath3D,[],2));
        yzPrjOverlay = cat(3,max(im,bw),im-bw,im-bw);
        imwrite(yzPrjOverlay, [outPath 'rPathOverlayYZ.tif'], 'tif', 'compression', 'none');
        
    end
    
end

end


function [pt1,pt2] = pickPointsFromProjectionImage(im)

% z-projection
imPrj = im2double(max(im,[],3));
% normalize the intensity range to [0,1]
imPrj = (imPrj - min(imPrj(:)))/(max(imPrj(:)) - min(imPrj(:)));

% manually select two tip points
imDisp=cat(3,imPrj,imPrj,imPrj);
figure(1),imshow(imDisp);
imageHandle = get(gca,'Children');

nPts = 0;
coords = zeros(2);
cent = zeros(size(imPrj));
while nPts < 2
    [xi,yi] = ginput(1);
    nPts = nPts + 1;
    coords(nPts,:) = round([xi yi]);
    cent(round(yi),round(xi))=1;
    centDisp=imdilate(cent,strel('disk',2));
    imDisp=cat(3,imPrj-centDisp,max(imPrj,centDisp),max(imPrj,centDisp));
    imDisp(imDisp < 0) = 0;
    set(imageHandle ,'CData',imDisp);
    drawnow;
end
pt1 = [coords(1,2),coords(1,1)];
pt2 = [coords(2,2),coords(2,1)];

end


function [pt1,pt2] = pickPointsFromColorImage(im)

% manually select two tip points
figure,imshow(im);
imageHandle = get(gca,'Children');

nPts = 0;
coords = zeros(2);
cent = zeros(size(im,1),size(im,2));
while nPts < 2
    [xi,yi] = ginput(1);
    nPts = nPts + 1;
    coords(nPts,:) = round([xi yi]);
    cent(round(yi),round(xi))=1;
    centDisp=imdilate(cent,strel('disk',2));
    imDisp=cat(3,max(im(:,:,1),centDisp),max(im(:,:,2),centDisp),max(im(:,:,3),centDisp));
    set(imageHandle ,'CData',imDisp);
    drawnow;
end
pt1 = [coords(1,2),coords(1,1)];
pt2 = [coords(2,2),coords(2,1)];

end


function [pt1,pt2] = findTipPointsOfAvgTrace(bwPaths3D,pt1s,pt2s)
% this function automatically finds the two tip points of the average trace
% of the three colored traces

%% 1st find the pair of tip points from the 3 colored traces on the opposite end that are the closest
% calculate the average position of all the pt1s
avgPt1 = mean(cat(1,pt1s{1},pt1s{2},pt1s{3}),1);
avgPt2 = mean(cat(1,pt2s{1},pt2s{2},pt2s{3}),1);

% find the pt1 that's closes to avgPt2
minDist = 1000;
for i = 1:3
    pt2i = pt2s{i};
    distVal = sqrt((avgPt1(1)-pt2i(1))^2 + (avgPt1(2)-pt2i(2))^2 + (avgPt1(3)-pt2i(3))^2);
    if (distVal < minDist)
        minDist = distVal;
        minInd1 = i;
    end
end
minPt1 = pt1s{minInd1};

% find the pt2 that's closes to avgPt1
minDist = 1000;
for i = 1:3
    pt1i = pt1s{i};
    distVal = sqrt((avgPt2(1)-pt1i(1))^2 + (avgPt2(2)-pt1i(2))^2 + (avgPt2(3)-pt1i(3))^2);
    if (distVal < minDist)
        minDist = distVal;
        minInd2 = i;
    end
end
minPt2 = pt2s{minInd2};

%% find the points on the other two traces that are closest to the minPts
% for pt1
group1 = zeros(3,3);
distmap = false(size(bwPaths3D{1}));
distmap(minPt1(1),minPt1(2),minPt1(3)) = true;
distmap = bwdist(distmap);
for i = 1:3
    if i == minInd1
        group1(i,:) = minPt1;
    else
        bwPath3D = bwPaths3D{i};
        tmp = distmap;
        tmp(~bwPath3D) = 1000;
        [~,ind] = min(tmp(:));
        [ymin,xmin,zmin] = ind2sub(size(bwPath3D),ind);
        group1(i,:) = [ymin,xmin,zmin];
    end
end

% for pt2
group2 = zeros(3,3);
distmap = false(size(bwPaths3D{1}));
distmap(minPt2(1),minPt2(2),minPt2(3)) = true;
distmap = bwdist(distmap);
for i = 1:3
    if i == minInd2
        group2(i,:) = minPt2;
    else
        bwPath3D = bwPaths3D{i};
        tmp = distmap;
        tmp(~bwPath3D) = 1000;
        [~,ind] = min(tmp(:));
        [ymin,xmin,zmin] = ind2sub(size(bwPath3D),ind);
        group2(i,:) = [ymin,xmin,zmin];
    end
end

%% take the average position of each group
pt1 = round(mean(group1,1));
pt2 = round(mean(group2,1));

end


function [blue,green,red] = loadImagesNeuralProcessesSpacialRelationship_TIF(filename)

r = bfopen(filename);
r = r{1};
blue = [];
green = [];
red = [];
for i = 1:size(r,1)
    tmp = r{i,2};
    nTmp = strfind(tmp,' C?=');
    cInd = str2num(tmp(nTmp+4));
    if cInd == 1
        blue = cat(3,blue,r{i,1});
    elseif cInd == 2
        green = cat(3,green,r{i,1});
    elseif cInd == 3
        red = cat(3,red,r{i,1});
    end
end

end


