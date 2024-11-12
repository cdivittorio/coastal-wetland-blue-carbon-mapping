%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Courtney Di Vittoriooo
% July 6 2023
% The purpose of this script is to read in GEE classification results
% and perform accuracy assessment, using test data for each fold of
% cross-validation.

% modified for classification scheme #2, which combines forest and scrub
% class

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%load matlab training data file, which contains all points
%find lat and lon that from full training that are not present in cv
%dataset and record type and year
close all
clear all

%training data folder - ENTER DIRECTORY WITH TRAINING DATA AND IMAGES HERE
cd('D:\Jacob Louie\ws8\cv1')
load('example_training_Barataria.mat')

%% ENTER FOLDER WHERE GEE CLASSIFICATION RESULTS ARE (IF IN A DIFFERENT FOLDER)
% You exported these from GEE - here is an example
%cd('D:\Jacob Louie\ws8\cv1')
%% get lat and lon of test pixels for each fold of cross validation
latFull = trainingBAR.full.latitude;
lonFull = trainingBAR.full.longitude;
yrFull = year(trainingBAR.full.date);
clFull = trainingBAR.full.type;
%start with cv1 %%%%%% CHANGE WHEN DOING OTHER CROSS-VALIDATION SETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
latTrain = trainingBAR.cv1.latitude;
lonTrain = trainingBAR.cv1.longitude;


k=0;
for j = 1:length(latFull)
   tmpLat = latFull(j);
   tmpLon = lonFull(j);
   tmpInd = find(latTrain == tmpLat & lonTrain == tmpLon);
   if isempty(tmpInd) == 1 %then part of test data
       k=k+1;
       latTest(k,1) = tmpLat;
       lonTest(k,1) = tmpLon;
       testYr(k,1) = yrFull(j);
       testCl(k,1) = clFull(j);
   end    
end
%make table
testInfo = table(lonTest,latTest,testYr,testCl,'VariableNames',...
    {'longitude','latitude','year','class'});
%% build table of lat and lon extents for each grid - same for all cv folds
% MODIFY FOR GRID ROWS AND COLUMNS IN YOUR WATERSHED
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Barataria shown as an example
gridNums = {'11','12','13','14','15',...
    '21','22','23','24','25','26','27',...
    '31','32','33','34','35','36','37','38',...
    '41','42','43','44','45','46','47','48','49',...
    '51','52','53','54','55','56','57','58','59'};
%loop that reads in files and grabs extents of grid - CHANGE FILE NAME FOR
%EACH CROSS-VALIDATION SET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:numel(gridNums)
   tmpFile = ['classificationProbs-cv1g',gridNums{j}]; 
   tmpInfo = geotiffinfo(tmpFile);
   lonMin(j,1) = tmpInfo.BoundingBox(1,1);
   lonMax(j,1) = tmpInfo.BoundingBox(2,1);
   latMin(j,1) = tmpInfo.BoundingBox(1,2);
   latMax(j,1) = tmpInfo.BoundingBox(2,2);
end
%make table of grid extents
gridExtents = table(gridNums',lonMin,lonMax,latMin,latMax,'VariableNames',...
    {'gridNum','lonMin','lonMax','latMin','latMax'});


%%  get resulting classification for each test point

%allocate memory
%geeClass(1:size(testInfo),1) = NaN;
nullVals = zeros(size(testInfo,1),1);
for j = 1:size(testInfo,1)
    tic
    %get lat and lon, year and class
    tmpLat = testInfo.latitude(j);
    tmpLon = testInfo.longitude(j);
    tmpYr = testInfo.year(j);
    tmpCl = testInfo.class(j);
    %find grid that pixels lies within
    tmpGridInd = find(tmpLon > gridExtents.lonMin & tmpLon < gridExtents.lonMax ...
        & tmpLat > gridExtents.latMin & tmpLat < gridExtents.latMax);
    tmpGrid = gridExtents.gridNum{tmpGridInd};
    %read in file - CHANGE FILE NAME FOR EACH CROSS VALIDATION DATASET
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmpFile = ['classificationProbs-cv1g',tmpGrid]; 
    tmpImg = geotiffread(tmpFile);
    tmpInfo = geotiffinfo(tmpFile);
    %create lat and lon grid
    height = tmpInfo.Height; % Integer indicating the height of the image in pixels
    width = tmpInfo.Width; % Integer indicating the width of the image in pixels
    [cols,rows] = meshgrid(1:width,1:height);
    [lonGrid,latGrid] = pix2map(tmpInfo.RefMatrix, rows, cols);
    %use grid locator function to find row and column wihtin image that
    %aligns with test pixel
    [tmpRow, tmpCol] = gridloc(tmpLat,tmpLon,latGrid,lonGrid);
    %get class probabilities for test year
    pixData(:,1) = tmpImg(tmpRow,tmpCol,:);
    %based on year of training, find position of image in stack
    layerInd = (6*(tmpYr-1996)+2:6*(tmpYr-1996)+6)';
    classProbs = pixData(layerInd);
    %if all values are null, data is missing and should skip
    if sum(isnan(classProbs)) == 5
        geeClass(j,1) = NaN;
        nullVals(j)=1;
        display(['image with null values for grid ',tmpGrid])
    else
        %most probable class
        geeClass(j,1) = find(classProbs == max(classProbs));
    end
    toc
end

%%

%switch gee class 4 to 5 and 5 to 6 so classes align in confusion matrix
geeClass(geeClass == 5) = 6;
geeClass(geeClass == 4) = 5;
%add GEE classification to table
%testInfo = addvars(testInfo,geeClass);
testInfo.geeClass = geeClass;
%% create confusion matrix 
%ccap in column, gee in row
classNums = [1 2 3 5 6];
for c = 1:5
    tmpInd = find(testInfo.class == classNums(c));
    tmpComp = testInfo.geeClass(tmpInd);
    % pixels that agree
    confMat(c,c) = length(tmpComp(tmpComp==classNums(c)));
    %false positive - gee classified pixels to this class (c) that should
    %not be classified as such
    tmpInd2 = find(testInfo.geeClass == classNums(c) & testInfo.class ~= classNums(c));
    %get values from test data that do not match 
    fpVals = testInfo.class(tmpInd2);
    for d = c+1:5
       %confMat(c,d) = 100*length(fpVals(fpVals == d))/length(find(testInfo.geeClass == c)); 
       confMat(c,d) = length(fpVals(fpVals == classNums(d))); 
    end
   %false negative - for test class pixels that were incorrectly labeled by gee, what did gee classify as instead?  
   for d = c+1:5
       %confMat(d,c) = 100*length(tmpComp(tmpComp == d))/length(tmpInd);
       confMat(d,c) = length(tmpComp(tmpComp == classNums(d)));
   end
end

figure
subplot(1,2,1)
confMat2 = confusionchart(testInfo.class,testInfo.geeClass,'RowSummary','row-normalized','ColumnSummary','column-normalized')

%% make class 6 = class 5 to make combined "others class
geeClass2 = testInfo.geeClass;
geeClass2(geeClass2 == 6) = 5;
trueClass2 = testInfo.class; 
trueClass2(trueClass2 == 6) = 5;

classNums = [1 2 3 5];
for c = 1:5
    tmpInd = find(trueClass2 == classNums(c));
    tmpComp = geeClass2(tmpInd);
    % pixels that agree
    confMat3(c,c) = length(tmpComp(tmpComp==classNums(c)));
    %false positive - gee classified pixels to this class (c) that should
    %not be classified as such
    tmpInd2 = find(geeClass2 == classNums(c) & trueClass2 ~= classNums(c));
    %get values from test data that do not match 
    fpVals = trueClass2(tmpInd2);
    for d = c+1:5
       %confMat(c,d) = 100*length(fpVals(fpVals == d))/length(find(testInfo.geeClass == c)); 
       confMat3(c,d) = length(fpVals(fpVals == classNums(d))); 
    end
   %false negative - for test class pixels that were incorrectly labeled by gee, what did gee classify as instead?  
   for d = c+1:5
       %confMat(d,c) = 100*length(tmpComp(tmpComp == d))/length(tmpInd);
       confMat3(d,c) = length(tmpComp(tmpComp == classNums(d)));
   end
end

%overall accuracy
100*sum(diag(confMat3))/sum(confMat3(:))

%use built-in matlab confusion matrix function
subplot(1,2,2)
confMat4 = confusionchart(trueClass2,geeClass2,'RowSummary','row-normalized','ColumnSummary','column-normalized')


%% SAVE RESULTS
%take screen shot of confusion matrices and add to sheet
% add overall accuracy and kappa index to summary table
%%%%% SAVE testInfo and confMat AS ws##cs#cv#.MAT FILE - RENAME FOR EACH
%%%%% CROSS-VALIDATION AND WATERSHED NUMBER 


%% Grid locator function

%gridloc - takes in file name, indices for band data to grab, 
% list of lat and lon coordinates, and gridsize and
%returns list of upper left rows and columns

function [rows,cols] = gridloc(samplat,samplon, latgrid, longrid)
    lat = latgrid;
    lon = longrid;
    for s = 1:length(samplat)
        latdiff = lat-samplat(s);
        londiff = lon-samplon(s);
        %min diff between lat and lon
        totdiff = abs(latdiff)+abs(londiff);
        [~, minloc] = min(abs(totdiff(:)));
        %get row and col info for location of min value - center
        [rowc(s,1),colc(s,1)]=ind2sub(size(totdiff),minloc);
    end
    rows = rowc;
    cols = colc;
end


