%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Courtney DiVittorio
% Aug 2023
% Now that the procedure has been established, it can be applied to the
% full images
% Steps
% 1 - create stacked annual classification maps that align with grid
% extents - for each grid find rows and columns that match annual map.
% create matlab data for each grid that can be pulled in concurrent with
% ccdc results
% 2 - iterate over ccdc results, 1 grid at a time
% 3 - for each pixel, grab ccdc probabilities and annual probabilities over time series 
% 4 - identify segment breaks and perform mixed class assignment and
% smoothing process
% Modified Dec 18 for more generalized use
%% 
clear all
% ws 51
gridNums = {'11','12','14','15','16',...
    '21','22','23','24','25','26',...
    '31','32','33','34','35','36',...
    '41','42','43','44','45','46',...
    '51','52','53','54','55','56',...]
    '61','62','63','64','65','66',...
    '71','72','73','74',...
    '81','82','83','84','85',...
    '93','94','95'};
%get lat and lon grid of annual maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE FOLDER FOR EACH WATERSHED 
cd('G:\My Drive\Research\NASA_OBB_Coastal_Marshes\Coastwide_Classification\Jacob_results\ws51_new\Annual')
info = geotiffinfo('classProbsMedMinMax-full2000.tif');
% get lat and lon grid
height = info.Height; % Integer indicating the height of the image in pixels
width = info.Width; % Integer indicating the width of the image in pixels
[cols,rows] = meshgrid(1:width,1:height);
[lonGridFull,latGridFull] = pix2map(info.RefMatrix, rows, cols);
%%
for gr = 1:length(gridNums)
    tic
    %go to ccdc results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE FOLDER FOR EACH WATERSHED 
    cd('G:\My Drive\Research\NASA_OBB_Coastal_Marshes\Coastwide_Classification\Jacob_results\ws51_new\CCDC')
    %import image
    [ccdcImg,R] = geotiffread(['classificationProbs-g',gridNums{gr},'.tif']);
    info = geotiffinfo(['classificationProbs-g',gridNums{gr},'.tif']);
    % get lat and lon grid
    height = info.Height; % Integer indicating the height of the image in pixels
    width = info.Width; % Integer indicating the width of the image in pixels
    [cols,rows] = meshgrid(1:width,1:height);
    [lonGrid,latGrid] = pix2map(info.RefMatrix, rows, cols);

    if latGrid(1,1) < latGrid(end,1)
        ccdcImg = flipud(ccdcImg);
        latGrid = flipud(latGrid);
    end
    %find row and col index that aligns with annual maps
    %upper left
    tmpDiff = abs(latGridFull-latGrid(1,1))+abs(lonGridFull-lonGrid(1,1));
    [minVal,minInd] = min(tmpDiff(:));
    [ulrow,ulcol] = ind2sub(size(tmpDiff),minInd);
    %upper right
    tmpDiff = abs(latGridFull-latGrid(1,end))+abs(lonGridFull-lonGrid(1,end));
    [minVal,minInd] = min(tmpDiff(:));
    [urrow,urcol] = ind2sub(size(tmpDiff),minInd);
     %lower left
    tmpDiff = abs(latGridFull-latGrid(end,1))+abs(lonGridFull-lonGrid(end,1));
    [minVal,minInd] = min(tmpDiff(:));
    [llrow,llcol] = ind2sub(size(tmpDiff),minInd);
    %lower right
    tmpDiff = abs(latGridFull-latGrid(end,end))+abs(lonGridFull-lonGrid(end,end));
    [minVal,minInd] = min(tmpDiff(:));
    [lrrow,lrcol] = ind2sub(size(tmpDiff),minInd);
    % annual maps full
    anImg = NaN(lrrow-ulrow+1,lrcol-ulcol+1,5*(2023-1984));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE FOLDER FOR EACH WATERSHED 
    cd('G:\My Drive\Research\NASA_OBB_Coastal_Marshes\Coastwide_Classification\Jacob_results\ws51_new\Annual')
    for yr = 1985:2023
        tmpImg = geotiffread(['classProbsMedMinMax-full',num2str(yr),'.tif']);
        anImg(:,:,(yr-1985)*5+1:(yr-1984)*5) = tmpImg(ulrow:lrrow,ulcol:lrcol,:);
    end
    %save to geotiff
    geotiffwrite(['classProbsMedMinMax-g',gridNums{gr},'.tif'],anImg,R)
    toc
end

disp('done with grid subsetting')
%% Use segment-based probabilities to classify
% for each grid, go through each pixel and grab ccdc and anual classifications probs
% perform classificaiton & assign mixed classes if below threshold
%
%class labels
classLabs = [1 2 3 4 5];
%mixed thresholds
mixedThresh = [0.47 0.40 0.41 0.32 0.34];

for gr = 1:length(gridNums)
    tic
    % read in annual classification probabilities for each grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE FOLDER FOR EACH WATERSHED
    cd('G:\My Drive\Research\NASA_OBB_Coastal_Marshes\Coastwide_Classification\Jacob_results\ws51_new\Annual')  
    anImg = geotiffread(['classProbsMedMinMax-g',gridNums{gr},'.tif']);
    %read in ccdc classification probabilities for each grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE FOLDER FOR EACH WATERSHED 
    cd('G:\My Drive\Research\NASA_OBB_Coastal_Marshes\Coastwide_Classification\Jacob_results\ws51_new\CCDC') 
    [ccdcImg,R] = geotiffread(['classificationProbs-g',gridNums{gr},'.tif']);
    %make empty classified image for memory allocation
    classImg = NaN(size(ccdcImg,1),size(ccdcImg,2),39);
    %change folder for storing results - might need to create this folder
    %in advance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE FOLDER FOR EACH WATERSHED 
    cd('G:\My Drive\Research\NASA_OBB_Coastal_Marshes\Coastwide_Classification\Jacob_results\ws51_new\matlab_output\classified_mixed') 
    %%%%%%%%%% NO CHANGES NEEDED BEYOND THIS POINT FOR THIS SECTION OF CODE
    %for each row and col
    for r = 1:size(ccdcImg,1)
        for c = 1:size(ccdcImg,2)
            clear ccdcPix 
            ccdcPix(:,1) = ccdcImg(r,c,:); %234 (39*6) - 1st layer is highest probability, followed by prob for each class
            anPix(:,1) = anImg(r,c,:); %195 (39*5) - does not contain 1st layer, just prob for each class
            %empty matrices for memory
            ccdcProbs = zeros(39,5);
            anProbs = zeros(39,5);
            anClass = zeros(39,1);
            ccdcClass = zeros(39,1);
            %rearrange so years in each row and probs in each column
            for yr = 1:39
                ccdcProbs(yr,:) = ccdcPix((yr-1)*6+2:yr*6); 
                anProbs(yr,:) = anPix((yr-1)*5+1:yr*5);
                %most likely class based on ccdc
                if sum(isnan(ccdcProbs(yr,:))) == 0
                    tmpMax = find(ccdcProbs(yr,:) == max(ccdcProbs(yr,:))); 
                    ccdcClass(yr,1) = tmpMax(1); %position not class
                else
                    ccdcClass(yr,1) = NaN;
                end
                if sum(isnan(anProbs(yr,:))) == 0
                    tmpMax = find(anProbs(yr,:) == max(anProbs(yr,:))); %might be more than 1
                    anClass(yr,1) = tmpMax(1); %position not class. take 1st class if more than 1
                else
                    anClass(yr,1) = NaN;
                end
            end
            %fill missing values
            ccdcClass = fillmissing(ccdcClass,'previous');
            ccdcClass = fillmissing(ccdcClass,'next');
            anClass = fillmissing(anClass,'previous');
            anClass = fillmissing(anClass,'next');
            %set ccdc as class in image
            classImg(r,c,:) = ccdcClass;
            for k = 1:5
                ccdcProbs(:,k) = fillmissing(ccdcProbs(:,k),'previous');
                ccdcProbs(:,k) = fillmissing(ccdcProbs(:,k),'next');
                anProbs(:,k) = fillmissing(anProbs(:,k),'previous');
                anProbs(:,k) = fillmissing(anProbs(:,k),'next');
            end
            %if below threshold relabel as mixed
            %find segment breaks
            locBreak = find(abs(diff(ccdcProbs(:,1)))>0); %represents end of segment but last segment not included
            if isempty(locBreak)==1
                %only 1 segment
                numSeg = 1;
                locBreak = size(ccdcProbs,1); %end of segment = 39
            else
                numSeg = length(locBreak)+1;
                locBreak(end+1) = size(ccdcProbs,1);
            end
            %loop over each segment to see if should be relabeled as mixed
            for k = 1:numSeg
                %grab probability of each class from ccdc
                tmpProbs = ccdcProbs(locBreak(k),:); %just look at 1st year b/c years that follow are same
                %find dominant class and probability
                domClass = find(tmpProbs == max(tmpProbs));
                domProb = max(tmpProbs);
                clear anSegProbs
                %test if less than threshold
                if isnan(domProb) == 1 %some pixels do not have data from ccdc 
                    %these should be corrected with spatial filter at end
                    classImg(r,c,s1:s2) = NaN;
                elseif domProb < mixedThresh(domClass(1)) %evaluate 
                    %get annual results over that segment
                    if k == 1
                        s1 = 1; s2 = locBreak(k);
                        anSeg = anClass(s1:s2,1);
                        anSegProbs(:,:) = anProbs(s1:s2,:);
                    else
                        s1 = locBreak(k-1)+1; s2 = locBreak(k);
                        anSeg = anClass(s1:s2,1);
                        anSegProbs(:,:) = anProbs(s1:s2,:);
                    end
                    %test if anSegProbs has single class
                    anSegClasses = unique(anSeg); %puts in order
                    if length(anSegClasses) == 1
                        %does class agree with CCDC? 
                        if anSeg(1) == domClass(1) %full class assigned
                            classImg(r,c,s1:s2) = classLabs(domClass(1)); %class instead of position
                        elseif anSeg(1) ~= domClass(1) 
                            % change to mixed
                            %use function that determines mixed class
                            classImg(r,c,s1:s2) = findMixed(anSeg(1),domClass(1));
                        end
                    elseif length(anSegClasses) == 2
                        %do classes match top 2 ccdc classes?
                        tmpProbs2 = tmpProbs;
                        tmpProbs2(domClass(1)) = 0;
                        domClass2 = find(tmpProbs2 == max(tmpProbs2));
                        if (anSegClasses(1) == domClass(1) | anSegClasses(1) == domClass2(1)) ...
                               && (anSegClasses(2) == domClass(1) | anSegClasses(2) == domClass2(1))
                            %they match
                            classImg(r,c,s1:s2) = findMixed(domClass(1),domClass2(1));
                        else %get average probability for each unique class and then take average of ccdc
                            for m = 1:size(anSegProbs,2)
                                tmpProbsAn(m) = mean(anSegProbs(:,m));
                            end
                            %get average
                            tmpProbsAvg = (tmpProbs+tmpProbsAn)./2; 
                            %use top two for mixed class
                            avgDom = find(tmpProbsAvg == max(tmpProbsAvg));
                            tmpProbsAvg2 = tmpProbsAvg;
                            tmpProbsAvg2(avgDom)=0;
                            avgDom2 = find(tmpProbsAvg2 == max(tmpProbsAvg2));
                            classImg(r,c,s1:s2) = findMixed(avgDom(1),avgDom2(1));
                        end
                    else %more than 2 - take average of both classificaitons
                        for m = 1:size(anSegProbs,2)
                            tmpProbsAn(m) = mean(anSegProbs(:,m));
                        end
                        %get average
                        tmpProbsAvg = (tmpProbs+tmpProbsAn)./2; 
                        %use top two for mixed class
                        avgDom = find(tmpProbsAvg == max(tmpProbsAvg));
                        tmpProbsAvg2 = tmpProbsAvg;
                        tmpProbsAvg2(avgDom)=0;
                        avgDom2 = find(tmpProbsAvg2 == max(tmpProbsAvg2));
                        classImg(r,c,s1:s2) = findMixed(avgDom(1),avgDom2(1));
                    end
                else %according to ccdc not mixed, but check annual results to see if they consistently disagree
                    %get annual results over that segment
                    clear anSegProbs
                    if k == 1
                        s1 = 1; s2 = locBreak(k);
                        anSeg = anClass(s1:s2,1);
                        anSegProbs(:,:) = anProbs(s1:s2,:);
                    else
                        s1 = locBreak(k-1)+1; s2 = locBreak(k);
                        anSeg = anClass(s1:s2,1);
                        anSegProbs(:,:) = anProbs(s1:s2,:);
                    end
                    %set values to full class for now and then consider if should
                    %be mixed
                    classImg(r,c,s1:s2) = classLabs(domClass(1));
                    %if segment is at least 3 years
                    if length(anSeg) > 2
                        %function that gets a value for how many steps the number stays at the same value
                        x = find(diff(anSeg)'); 
                        n = [x numel(x)] - [0 x];
                        y = arrayfun(@(X) X-1:-1:0, n , 'un',0);
                        consec = cat(2,y{:});
                        %note - help from mathworks online
                        %find values 2 or greater
                        tmpInd = find(consec>=2);
                        %for each segment check and see if different than ccdc
                        %class and change segment to mixed
                        for m = 1:length(tmpInd)
                            if anSeg(tmpInd(m)) ~= domClass(1)
                                s1a = s1+tmpInd(m)-1; %starting point of sub-segment
                                s2a = s1+tmpInd(m)+consec(tmpInd(m))-1; %end point of sub-segment
                                classImg(r,c,s1a:s2a) = findMixed(domClass(1),anSeg(tmpInd(m)));
                            end
                        end
                    end
                end
            end

        end
    end
    toc
    save(['classImg-g',gridNums{gr},'.mat'],'classImg');
end
disp('done with mixed classification')
%%  smoothing and change characterization

% create look-up table with acceptable transitions, where there is a common class between segments
% 11 - water/emergent
% 12 - water/forest-scrub
% 13 - water/other
% 14 - emergent/forest-scrub
% 15 - emergent/other
% 16 - forest-scrub/other
%

commonClassComb = [11 12; 11 13; 12 13; ...
    11 14; 11 15; 14 15; 12 14; 12 16; 14 16; ...
    13 15; 13 16; 15 16; 1 11; 1 12; 1 13;...
    2 11; 2 14; 2 15; 3 12; 3 14; 3 16;...
    4 13; 4 15; 4 16; 4 5; 5 13; 5 15; 5 16];
%switch columns and append
commonClassComb(end+1:size(commonClassComb,1)*2,1) = commonClassComb(:,2);
commonClassComb(size(commonClassComb,1)/2+1:end,2) = commonClassComb(1:size(commonClassComb,1)/2,1);

%create look-up table for mixed classes and the classes they contain, where
%column represents class - water/emerg/forest-scrub/other
mixedClassesTab = [1 2 3 4 5; 11 11 12 13 13; 12 14 14 15 15; 13 15 16 16 16];

% load one grid at a time
% look at full time series and force transitions to be logical
for gr = 1:length(gridNums)
    tic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE FOLDER FOR EACH WATERSHED 
    cd('G:\My Drive\Research\NASA_OBB_Coastal_Marshes\Coastwide_Classification\Jacob_results\ws51_new\matlab_output\classified_mixed') 
    load(['classImg-g',gridNums{gr},'.mat']);
    changeType = NaN(size(classImg,1),size(classImg,2));
    years = (1985:2023)';
    %get probabilities from annual and ccdc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE FOLDER FOR EACH WATERSHED
    cd('G:\My Drive\Research\NASA_OBB_Coastal_Marshes\Coastwide_Classification\Jacob_results\ws51_new\Annual')  
    anImg = geotiffread(['classProbsMedMinMax-g',gridNums{gr},'.tif']);
    %ccdc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE FOLDER FOR EACH WATERSHED
    cd('G:\My Drive\Research\NASA_OBB_Coastal_Marshes\Coastwide_Classification\Jacob_results\ws51_new\CCDC')  
    ccdcImg = geotiffread(['classificationProbs-g',gridNums{gr},'.tif']);
   
    %%%%%%%%%%%%%%%%%%%CHANGE FOLDER FOR EACH WATERSHED - MAKE THIS FOLDER
    cd('G:\My Drive\Research\NASA_OBB_Coastal_Marshes\Coastwide_Classification\Jacob_results\ws51_new\matlab_output\classified_smoothed')  
    %%%%%%%%%%%%%%%%%%%%%%%% NO MORE CHANGES NEEDED IN THIS SECTION OF CODE
    %new time series of smoothed classifications
    for r = 1:size(classImg,1)
        %tic
        for c = 1:size(classImg,2)
            %tic
            clear pixclasses
            pixclasses(1,:) = classImg(r,c,:);
            %change 5 to 4 so does not count changes between 4 and 5 as change
            pixclasses(pixclasses == 5) = 4;
            clear changeinfo changeClasses
            changeinfo(1) = 0; %cant have change on 1st yr
            changeinfo(2:length(pixclasses),:) = diff(pixclasses);
            changeinfo(changeinfo ~=0) = 1;
            %get years where changes occur and count total number of
            %changes 
            changeYrs = years(changeinfo==1);
            totChange = length(changeYrs);
            changeloc = find(changeinfo == 1);
            %find classes for each change
            if isempty(changeloc) == 0
                changeClasses(1) = pixclasses(changeloc(1)-1);
                changeClasses(2:length(changeloc)+1) = pixclasses(changeloc);
                mixedClasses = changeClasses(changeClasses>10);
                fullClasses = changeClasses(changeClasses<10);
                mixedChangeOrder = find(changeClasses>10);
            end
            %if no changes then give value of 0
            if sum(changeinfo) == 0
                changeType(r,c,1) = 0; %NO CHANGE
                %numPermChanges(j,1) = 0;
                %numTempChanges(j,1) = 0;
            elseif sum(changeinfo)>0 %at least 1 change 
                %is there a change from 1 full class to another full class? abrupt
                %change. If mixed class between then changes are gradual
                %could have both abrupt and gradual changes.
        
                % Are changes between mixedIs class at end the same as class at beginning? If so, then
                % temporary change
                %create temp/perm info
                changeloc = find(changeinfo == 1);
                %find classes for each change
                changeClasses(1) = pixclasses(changeloc(1)-1);
                changeClasses(2:length(changeloc)+1) = pixclasses(changeloc);
                mixedClasses = changeClasses(changeClasses>10);
                fullClasses = changeClasses(changeClasses<10);
        
                %if all classes greater than 10, then mixed change only
                if length(mixedClasses) == length(changeClasses)
                    %determine if 1st and last class the same
                    if changeClasses(1) == changeClasses(end)
                        %mixed change = temporary
                        changeType(r,c,1) = 1;
                    else
                        %mixed change - permanent
                        changeType(r,c,1) = 2;
                    end
                    %check to see if common class between mixed segments
                    %if not, find probabilities of each class over both segments
                    %class with highest overall probability should be common class
                    %for segment that needs to change, use other dominant class for
                    %mixed class
                    for k = 1:length(mixedClasses)-1
                        c1 = mixedClasses(k);
                        c2 = mixedClasses(k+1);
                        commonInd =  ismember([c1 c2], commonClassComb, 'rows');
                        if commonInd == 0 %segment needs to change
                            clear ccdcProbs anProbs
                            %use probabilities from image
                            
                            ccdcPix(:,1) = ccdcImg(r,c,:); %234 (39*6)
                            anPix(:,1) = anImg(r,c,:); %195 (39*5)
                            %function that organizes and fills missing
                            %values
                            [ccdcProbs, anProbs] = getProbs(ccdcPix,anPix);
                            %segment start and end
                            if length(mixedClasses) == 2
                                %look at full segment
                                %take average of rows
                                segProbs = mean((ccdcProbs+anProbs)./2,1);
                                s1Probs = mean((ccdcProbs(1:changeloc-1,:)+anProbs(1:changeloc-1,:))./2,1);
                                s2Probs = mean((ccdcProbs(changeloc:end,:)+anProbs(changeloc:end,:))./2,1);
                            elseif length(mixedClasses) > 2
                                if k == 1 %1st and 2nd segment
                                    s1a = 1; s1b = changeloc(k)-1;
                                    s2a = changeloc(k); s2b = changeloc(k+1)-1;
                                elseif k > 1 & k == length(mixedClasses)-1 %end of segment
                                    s1a = changeloc(k-1); s1b = changeloc(k)-1;
                                    s2a = changeloc(k); s2b = 39;
                                else %middle segments
                                    s1a = changeloc(k-1); s1b = changeloc(k)-1;
                                    s2a = changeloc(k); s2b = changeloc(k+1)-1;
                                end
                                s1Probs = mean((ccdcProbs(s1a:s1b,:)+anProbs(s1a:s1b,:))./2,1);
                                s2Probs = mean((ccdcProbs(s2a:s2b,:)+anProbs(s2a:s2b,:))./2,1);
                                segProbs = mean((ccdcProbs(s1a:s2b,:)+anProbs(s1a:s2b,:))./2,1);
                            end
                            %find max probability 
                            domClassFull = find(segProbs == max(segProbs));
                            domClass1 = find(s1Probs == max(s1Probs));
                            domClass2 = find(s2Probs == max(s2Probs));
                            %if segment 1 does not have one of dominant full
                            %classes then switch
                            if ismember(mixedClasses(k),mixedClassesTab(:,domClassFull(1))) == 0
                                if domClassFull(1) == domClass1(1)
                                    s1Probs(domClass1(1)) = 0;
                                    domClass1b = find(s1Probs == max(s1Probs));
                                    classImg(r,c,s1a:s1b) = findMixed(domClassFull(1),domClass1b(1));
                                else
                                    classImg(r,c,s1a:s1b) = findMixed(domClassFull(1),domClass1(1));
                                end
                                %display(['change for row',num2str(r),' col ',num2str(c)])
                            end
                            %repeat for 2nd segment
                            if ismember(mixedClasses(k+1),mixedClassesTab(:,domClassFull(1))) == 0
                                if domClassFull(1) == domClass2(1)
                                    s2Probs(domClass2(1)) = 0;
                                    domClass2b = find(s2Probs == max(s2Probs));
                                    classImg(r,c,s2a:s2b) = findMixed(domClassFull(1),domClass2b(1));
                                else
                                    classImg(r,c,s2a:s2b) = findMixed(domClassFull(1),domClass2(1));
                                end
                                %display(['change for row',num2str(r),' col ',num2str(c)])
                            end
                        end
                    end
                end 
                
                %there are full classes present
        
                %is there an abrupt transition?
                %if so, is there also a gradual transition?
                %if so make sure that gradual transition aligns with full class labels.
                % and change is abrupt and gradual - is the change temporary or permanent?
                %if only abrupt transitions, then it is abrupt full change -
                %determine if temp/perm.
                clear ccdcProbs anProbs
                ccdcPix(:,1) = ccdcImg(r,c,:); %234 (39*6)
                anPix(:,1) = anImg(r,c,:); %195 (39*5)
                %function that organizes and fills missing
                %values
                [ccdcProbs, anProbs] = getProbs(ccdcPix,anPix);
                if length(fullClasses) > 1 & length(unique(fullClasses))>1 %more than 1 full class
                    clear locFull
                    fullClassesUnique = unique(fullClasses,'stable'); %do not sort
                    for k = 1:length(fullClassesUnique)
                        tmpLoc = find(changeClasses == fullClassesUnique(k)); %records where full classes are
                        locFull(k) = tmpLoc(1); %first occurence of class
                    end
                    diffFull = diff(locFull); %to see if next to each other
                    if length(find(diffFull == 1))>0 %then abrupt transition present
                        %is there also a gradual transition?
                        if length(mixedClasses) > 0
                            %gradual also present
                            %  APPLY SMOOTHING SO MIXED/FULL CLASS TRANSITIONS HAVE COMMON
                            %  CLASS
                            %need to find segment start for mixed class
                            mixedChangeOrder = find(changeClasses>10); %index for change class that is mixed
                            %mixedSegs = changeloc(find(changeClasses>10));
                            for k = 1:length(mixedClasses)
                                %find segment before or after
                                if mixedChangeOrder(k) == 1 %first segment is mixed so take following segment as comparison
                                    c1 = changeClasses(1); c2 = changeClasses(2);
                                    s1a = 1; s1b = changeloc(1)-1;
                                    if length(changeloc) == 1 %only 1 change
                                        s2a = changeloc(1); s2b = 39;
                                    else 
                                        s2a = changeloc(1); s2b = changeloc(2)-1;
                                    end
                                    %GET PROBS AND ASSIGN CLASS
                                    s1Probs = mean((ccdcProbs(s1a:s1b,:)+anProbs(s1a:s1b,:))./2,1);
                                    s2Probs = mean((ccdcProbs(s2a:s2b,:)+anProbs(s2a:s2b,:))./2,1);
                                    segProbs = mean((ccdcProbs(s1a:s2b,:)+anProbs(s1a:s2b,:))./2,1);
                                    domClassFull = find(segProbs == max(segProbs));
                                    domClass1 = find(s1Probs == max(s1Probs));
                                    domClass2 = find(s2Probs == max(s2Probs));
                                    %if segment 1 does not have one of dominant full
                                    %classes then switch 
                                    if ismember(c1,mixedClassesTab(:,domClassFull(1))) == 0 & c2 > 10 
                                        if domClassFull(1) == domClass1(1)
                                            s1Probs(domClass1(1)) = 0;
                                            domClass1b = find(s1Probs == max(s1Probs));
                                            classImg(r,c,s1a:s1b) = findMixed(domClassFull(1),domClass1b(1));
                                        else
                                            classImg(r,c,s1a:s1b) = findMixed(domClassFull(1),domClass1(1));
                                        end
                                        %display(['change for row ',num2str(r), 'col ',num2str(c)])
                                    elseif ismember(c1,mixedClassesTab(:,domClassFull(1))) == 0 & c2 < 10
                                        if domClass1(1) == c2
                                            s1Probs(domClass1(1)) = 0;
                                            domClass1b = find(s1Probs == max(s1Probs));
                                            classImg(r,c,s1a:s1b) = findMixed(c2,domClass1b(1));
                                        else
                                            classImg(r,c,s1a:s1b) = findMixed(c2,domClass1(1));
                                        end
                                         %display(['change for row ',num2str(r), 'col ',num2str(c)])
                                    end
                                elseif mixedChangeOrder(k) == length(changeClasses) %last segment is mixed
                                    c1 = changeClasses(mixedChangeOrder(k)-1); 
                                    c2 = changeClasses(mixedChangeOrder(k));
                                    %if segment length = 1 for last 2
                                    %segments
                                    if changeloc(end-1)+1 == 39
                                        s1a = 38; s1b = 38; 
                                        s2a = 39; s2b = 39;
                                    %if first segment has length of 1
                                    elseif changeloc(end)-changeloc(end-1) == 1
                                        s1a = changeloc(end-1);
                                        s1b = changeloc(end-1); 
                                        s2a = changeloc(end);
                                        s2b = 39;
                                    else
                                        s1a = changeloc(end-1)+1;
                                        s1b = changeloc(end)-1; %end of last segment
                                        s2a = changeloc(end);
                                        s2b = 39;
                                    end
                                    %GET PROBS AND ASSIGN CLASS - COULD MAKE A
                                    %FUNCTION
                                    s1Probs = mean((ccdcProbs(s1a:s1b,:)+anProbs(s1a:s1b,:))./2,1);
                                    s2Probs = mean((ccdcProbs(s2a:s2b,:)+anProbs(s2a:s2b,:))./2,1);
                                    segProbs = mean((ccdcProbs(s1a:s2b,:)+anProbs(s1a:s2b,:))./2,1);
                                    domClassFull = find(segProbs == max(segProbs));
                                    domClass1 = find(s1Probs == max(s1Probs));
                                    domClass2 = find(s2Probs == max(s2Probs));
                                    %if segment 1 does not have one of dominant full
                                    %classes then switch
                                    if ismember(mixedClasses(k),mixedClassesTab(:,domClassFull(1))) == 0 & c1 > 10 %previous is also mixed
                                        if domClassFull(1) == domClass1(1) % if both have same dominant class find next most probably class
                                            s1Probs(domClass1(1)) = 0;
                                            domClass1b = find(s1Probs == max(s1Probs));
                                            classImg(r,c,s1a:s1b) = findMixed(domClassFull(1),domClass1b(1));
                                        else
                                            classImg(r,c,s1a:s1b) = findMixed(domClassFull(1),domClass1(1));
                                        end
                                    end
                                    if c1 < 10 
                                        if ismember(mixedClasses(k),mixedClassesTab(:,c1)) == 0 
                                            if c1 == domClass2(1)
                                                s2Probs(domClass2(1)) = 0; 
                                                domClass2b = find(s2Probs == max(s2Probs));
                                                classImg(r,c,s2a:s2b) = findMixed(domClass2b(1),c1);
                                            else
                                                classImg(r,c,s2a:s2b) = findMixed(domClass2(1),c1);
                                            end
                                        end
                                    end
                                else %middle segment,
                                    c1 = changeClasses(mixedChangeOrder(k)-1); 
                                    c2 = changeClasses(mixedChangeOrder(k));
                                    c3 = changeClasses(mixedChangeOrder(k)+1);
                                    s1a = changeloc(mixedChangeOrder(k)-1);
                                    s1b = changeloc(mixedChangeOrder(k))-1;
                                    s2a = changeloc(mixedChangeOrder(k));
                                    if length(changeloc) == mixedChangeOrder(k) %next segment is last segment
                                        s2b = 39;
                                    else
                                        s2b = changeloc(mixedChangeOrder(k)+1)-1;
                                    end
                                    %if c1 and c3 are full class, then c2 should be
                                    %mixed between them, but if full classes are
                                    %the same, then just make sure it has perm
                                    %class and 2nd most likely class for segment
                                    %
                                    if c1 < 10 &  c3 < 10 & c1 ~= c3
                                        classImg(r,c,s2a:s2b) = findMixed(c1,c3);
                                    elseif c1 < 10 & c1 == c3
                                        s2Probs = mean((ccdcProbs(s2a:s2b,:)+anProbs(s2a:s2b,:))./2,1);
                                        domClass2a = find(s2Probs == max(s2Probs));
                                        s2Probs(domClass2a(1)) = 0;
                                        domClass2b = find(s2Probs == max(s2Probs));
                                        if c1 ~= domClass2a(1)
                                             classImg(r,c,s2a:s2b) = findMixed(c1,domClass2a(1));
                                        else
                                             classImg(r,c,s2a:s2b) = findMixed(c1,domClass2b(1));
                                        end
                                    elseif c1 < 10 & c3 > 10 %following class is mixed and will be corrected later, so just look at previous
                                        s2Probs = mean((ccdcProbs(s2a:s2b,:)+anProbs(s2a:s2b,:))./2,1);
                                        domClass2 = find(s2Probs == max(s2Probs));
                                        if domClass2(1) == c1
                                            s2Probs(domClass2(1)) = 0;
                                            domClass2b = find(s2Probs == max(s2Probs));
                                            classImg(r,c,s2a:s2b) = findMixed(c1,domClass2b(1));
                                        else
                                            classImg(r,c,s2a:s2b) = findMixed(c1,domClass2(1));
                                        end
                                    else %previous and following mixed, align with previous
                                        s1Probs = mean((ccdcProbs(s1a:s1b,:)+anProbs(s1a:s1b,:))./2,1);
                                        s2Probs = mean((ccdcProbs(s2a:s2b,:)+anProbs(s2a:s2b,:))./2,1);
                                        segProbs = mean((ccdcProbs(s1a:s2b,:)+anProbs(s1a:s2b,:))./2,1);
                                        domClassFull = find(segProbs == max(segProbs));
                                        domClass1 = find(s1Probs == max(s1Probs));
                                        domClass2 = find(s2Probs == max(s2Probs));
                                        if ismember(mixedClasses(k),mixedClassesTab(:,domClassFull(1))) == 0 
                                            if domClassFull(1) == domClass2(1)
                                                s2Probs(domClass2)=0;
                                                domClass2b = find(s2Probs == max(s2Probs));
                                                classImg(r,c,s1a:s1b) = findMixed(domClass2b(1),domClassFull(1));
                                            else
                                                classImg(r,c,s1a:s1b) = findMixed(domClass2(1),domClassFull(1));
                                            end
                                            
                                        end
                                    end
        
                                end
        
                            end
                            % ASSIGN CHANGE TYPE Temporary or permanent?
                            if changeClasses(1) == changeClasses(end)
                                %abrupt and gradual change - temporary
                                changeType(r,c,1) = 9;
                            else
                                %abrupt and gradual change - permanent
                                changeType(r,c,1) = 10;
                            end
                        else %gradual not present, abrupt only
                             % Temporary or permanent?
                            if changeClasses(1) == changeClasses(end)
                                %abrupt change - temporary
                                changeType(r,c,1) = 7;
                            else
                                %abrupt change - permanent
                                changeType(r,c,1) = 8;
                            end
                        end
                    else % more than 1 full class but not next to each other so mixed class between
                        %gradual full change - determine temp/perm. 
                        %  APPLY SMOOTHING -
                        mixedChangeOrder = find(changeClasses>10); %index for change class that is mixed
                        %mixedSegs = changeloc(find(changeClasses>10));
                        for k = 1:length(mixedClasses)
                            %find segment before or after
                            if mixedChangeOrder(k) == 1 %first segment is mixed so take following segment as comparison & change 1st
                                c1 = changeClasses(1); c2 = changeClasses(2);
                                s1a = 1; s1b = changeloc(1)-1;
                                if length(changeloc) == 1 %only 1 change
                                    s2a = changeloc(1); s2b = 39;
                                else 
                                    s2a = changeloc(1); s2b = changeloc(2)-1;
                                end
                                %GET PROBS AND ASSIGN CLASS
                                s1Probs = mean((ccdcProbs(s1a:s1b,:)+anProbs(s1a:s1b,:))./2,1);
                                s2Probs = mean((ccdcProbs(s2a:s2b,:)+anProbs(s2a:s2b,:))./2,1);
                                segProbs = mean((ccdcProbs(s1a:s2b,:)+anProbs(s1a:s2b,:))./2,1);
                                domClassFull = find(segProbs == max(segProbs));
                                domClass1 = find(s1Probs == max(s1Probs));
                                domClass2 = find(s2Probs == max(s2Probs));
                                %if segment 1 does not have one of dominant full
                                %classes then switch to class that does
                                if ismember(c1,mixedClassesTab(:,domClassFull(1))) == 0 & c2 > 10 
                                    if domClassFull(1) == domClass1(1)
                                        s1Probs(domClass1(1)) = 0;
                                        domClass1b = find(s1Probs == max(s1Probs));
                                        classImg(r,c,s1a:s1b) = findMixed(domClassFull(1),domClass1b(1));
                                    else
                                        classImg(r,c,s1a:s1b) = findMixed(domClassFull(1),domClass1(1));
                                    end
                                elseif ismember(c1,mixedClassesTab(:,domClassFull(1))) == 0 & c2 < 10
                                    if domClass1(1) == c2
                                        s1Probs(domClass1(1)) = 0;
                                        domClass1b = find(s1Probs == max(s1Probs));
                                        classImg(r,c,s1a:s1b) = findMixed(c2,domClass1b(1));
                                    else
                                        classImg(r,c,s1a:s1b) = findMixed(c2,domClass1(1));
                                    end
                                end
                            elseif mixedChangeOrder(k) == length(changeClasses) %last segment
                                c1 = changeClasses(mixedChangeOrder(k)-1); 
                                c2 = changeClasses(mixedChangeOrder(k));
                                if  length(changeloc) == 1 %only 1 change
                                    s1a = 1;
                                    s1b = changeloc(1)-1;
                                    s2a = changeloc(1);
                                    s2b = 39;
                                else
                                    s1a = changeloc(mixedChangeOrder(k)-2);
                                    s1b = changeloc(mixedChangeOrder(k)-1)-1;
                                    s2a = changeloc(mixedChangeOrder(k)-1);
                                    if mixedChangeOrder(k) > length(changeloc)
                                        s2b = 39;
                                    else
                                        s2b = changeloc(mixedChangeOrder(k))-1;
                                    end
                                end
                           
                                %GET PROBS AND ASSIGN CLASS - COULD MAKE A
                                %FUNCTION
                                s1Probs = mean((ccdcProbs(s1a:s1b,:)+anProbs(s1a:s1b,:))./2,1);
                                s2Probs = mean((ccdcProbs(s2a:s2b,:)+anProbs(s2a:s2b,:))./2,1);
                                segProbs = mean((ccdcProbs(s1a:s2b,:)+anProbs(s1a:s2b,:))./2,1);
                                domClassFull = find(segProbs == max(segProbs));
                                domClass1 = find(s1Probs == max(s1Probs));
                                domClass2 = find(s2Probs == max(s2Probs));
                                %if segment 1 does not have one of dominant full
                                %classes then switch
                                if ismember(mixedClasses(k),mixedClassesTab(:,domClassFull(1))) == 0 & c1 > 10 %previous is also mixed
                                    if domClassFull(1) == domClass1(1)
                                        segProbs(domClass1(1)) = 0;
                                        domClass1b = find(segProbs == max(segProbs));
                                        classImg(r,c,s1a:s1b) = findMixed(domClassFull(1),domClass1b(1));
                                    else
                                        classImg(r,c,s1a:s1b) = findMixed(domClassFull(1),domClass1(1));
                                    end
                                end
                                if c1 < 10
                                    if ismember(mixedClasses(k),mixedClassesTab(:,domClassFull(1))) == 0 %prev is full class
                                        if domClass2(1) == c1
                                            s2Probs(domClass2(1)) = 0;
                                            domClass2b = find(s2Probs == max(s2Probs));
                                            classImg(r,c,s2a:s2b) = findMixed(domClass2b(1),c1);
                                        else
                                            classImg(r,c,s2a:s2b) = findMixed(domClass2(1),c1);
                                        end
                                    end
                                end
                            else %middle segment,
                                c1 = changeClasses(mixedChangeOrder(k)-1); 
                                c2 = changeClasses(mixedChangeOrder(k));
                                c3 = changeClasses(mixedChangeOrder(k)+1);
                                s1a = changeloc(mixedChangeOrder(k)-1);
                                s1b = changeloc(mixedChangeOrder(k))-1;
                                s2a = changeloc(mixedChangeOrder(k));
                                if length(changeloc) == mixedChangeOrder(k) %next segment is last segment
                                    s2b = 39;
                                else
                                    s2b = changeloc(mixedChangeOrder(k)+1)-1;
                                end
                                %if c1 and c3 are full class, then c2 should be
                                %mixed between them, but if full classes are
                                %the same, then just make sure it has perm
                                %class and 2nd most likely class for segment
                                %
                                if c1 < 10 &  c3 < 10 & c1 ~= c3
                                    classImg(r,c,s2a:s2b) = findMixed(c1,c3);
                                elseif c1 < 10 & c1 == c3
                                    s2Probs = mean((ccdcProbs(s2a:s2b,:)+anProbs(s2a:s2b,:))./2,1);
                                    domClass2a = find(s2Probs == max(s2Probs));
                                    s2Probs(domClass2a(1)) = 0;
                                    domClass2b = find(s2Probs == max(s2Probs));
                                    if c1 ~= domClass2a(1)
                                         classImg(r,c,s2a:s2b) = findMixed(c1,domClass2a(1));
                                    else
                                         classImg(r,c,s2a:s2b) = findMixed(c1,domClass2b(1)); %LEFT OFF HERE
                                    end
                                elseif c1 < 10 & c3 > 10 %following class is mixed and will be corrected later, so just look at previous
                                    s2Probs = mean((ccdcProbs(s2a:s2b,:)+anProbs(s2a:s2b,:))./2,1);
                                    domClass2 = find(s2Probs == max(s2Probs));
                                    if domClass2(1) == c1
                                        s2Probs(domClass2(1)) = 0;
                                        domClass2b = find(s2Probs == max(s2Probs));
                                        classImg(r,c,s2a:s2b) = findMixed(c1,domClass2b(1));
                                    else
                                        classImg(r,c,s2a:s2b) = findMixed(c1,domClass2(1));
                                    end
                                else %previous and following mixed, align with previous
                                    s1Probs = mean((ccdcProbs(s1a:s1b,:)+anProbs(s1a:s1b,:))./2,1);
                                    s2Probs = mean((ccdcProbs(s2a:s2b,:)+anProbs(s2a:s2b,:))./2,1);
                                    segProbs = mean((ccdcProbs(s1a:s2b,:)+anProbs(s1a:s2b,:))./2,1);
                                    domClassFull = find(segProbs == max(segProbs));
                                    domClass1 = find(s1Probs == max(s1Probs));
                                    domClass2 = find(s2Probs == max(s2Probs));
                                    if ismember(mixedClasses(k),mixedClassesTab(:,domClassFull(1))) == 0 
                                        if domClassFull(1) == domClass2(1)
                                            s2Probs(domClass2(1))=0;
                                            domClass2b = find(s2Probs == max(s2Probs));
                                            classImg(r,c,s1a:s1b) = findMixed(domClass2b(1),domClassFull(1));
                                        else
                                            classImg(r,c,s1a:s1b) = findMixed(domClass2(1),domClassFull(1));
                                        end
                                    end
                                end
        
                            end
                        end
                        if changeClasses(1) == changeClasses(end)
                            %gradual full change - temporary
                            changeType(r,c,1) = 5;
                        else
                            %gradual full change - permanent
                            changeType(r,c,1) = 6;
                        end
                    end
                else %only 1 full class and change is present
                    %so gradual mixed change 
                    %  APPLY SMOOTHING 
                    for k = 1:length(mixedClasses)
                        %find segment before or after
                        if mixedChangeOrder(k) == 1 %first segment is mixed so take following segment as comparison
                            c1 = changeClasses(1); c2 = changeClasses(2);
                            s1a = 1; s1b = changeloc(1)-1;
                            if length(changeloc) == 1 %only 1 change
                                s2a = changeloc(1); s2b = 39;
                            else 
                                s2a = changeloc(1); s2b = changeloc(2)-1;
                            end
                            %GET PROBS AND ASSIGN CLASS
                            s1Probs = mean((ccdcProbs(s1a:s1b,:)+anProbs(s1a:s1b,:))./2,1);
                            s2Probs = mean((ccdcProbs(s2a:s2b,:)+anProbs(s2a:s2b,:))./2,1);
                            segProbs = mean((ccdcProbs(s1a:s2b,:)+anProbs(s1a:s2b,:))./2,1);
                            domClassFull = find(segProbs == max(segProbs));
                            domClass1 = find(s1Probs == max(s1Probs));
                            domClass2 = find(s2Probs == max(s2Probs));
                            %if segment 1 does not have one of dominant full
                            %classes then switch to class that does
                            if ismember(c1,mixedClassesTab(:,domClassFull(1))) == 0 & c2 > 10 
                                if domClassFull(1) == domClass1(1)
                                    s1Probs(domClass1(1)) = 0;
                                    domClass1b = find(s1Probs == max(s1Probs));
                                    classImg(r,c,s1a:s1b) = findMixed(domClassFull(1),domClass1b(1));
                                else
                                    classImg(r,c,s1a:s1b) = findMixed(domClassFull(1),domClass1(1));
                                end
                            elseif ismember(c1,mixedClassesTab(:,domClassFull(1))) == 0 & c2 < 10
                                if domClass1(1) == c2
                                    s1Probs(domClass1(1)) = 0;
                                    domClass1b = find(s1Probs == max(s1Probs));
                                    classImg(r,c,s1a:s1b) = findMixed(c2,domClass1b(1));
                                else
                                    classImg(r,c,s1a:s1b) = findMixed(c2,domClass1(1));
                                end
                            end
                        elseif mixedChangeOrder(k) == length(changeClasses) %last segment
                            c1 = changeClasses(mixedChangeOrder(k)-1); 
                            c2 = changeClasses(mixedChangeOrder(k));
                            if  length(changeloc) == 1 %only 1 change
                                s1a = 1;
                                s1b = changeloc(1)-1;
                                s2a = changeloc(1);
                                s2b = 39;
                            else
                                s1a = changeloc(mixedChangeOrder(k)-2);
                                s1b = changeloc(mixedChangeOrder(k)-1)-1;
                                s2a = changeloc(mixedChangeOrder(k)-1);
                                if mixedChangeOrder(k) > length(changeloc)
                                    s2b = 39;
                                else
                                    s2b = changeloc(mixedChangeOrder(k))-1;
                                end
                            end
                            %GET PROBS AND ASSIGN CLASS - COULD MAKE A
                            %FUNCTION
                            s1Probs = mean((ccdcProbs(s1a:s1b,:)+anProbs(s1a:s1b,:))./2,1);
                            s2Probs = mean((ccdcProbs(s2a:s2b,:)+anProbs(s2a:s2b,:))./2,1);
                            segProbs = mean((ccdcProbs(s1a:s2b,:)+anProbs(s1a:s2b,:))./2,1);
                            domClassFull = find(segProbs == max(segProbs));
                            domClass1 = find(s1Probs == max(s1Probs));
                            domClass2 = find(s2Probs == max(s2Probs));
                            %if segment 1 does not have one of dominant full
                            %classes then switch
                            if ismember(mixedClasses(k),mixedClassesTab(:,domClassFull(1))) == 0 & c1 > 10 %previous is also mixed
                                if domClassFull(1) == domClass1(1)
                                    s1Probs(domClass1(1)) = 0;
                                    domClass1b = find(s1Probs == max(s1Probs));
                                    classImg(r,c,s1a:s1b) = findMixed(domClassFull(1),domClass1b(1));
                                else
                                    classImg(r,c,s1a:s1b) = findMixed(domClassFull(1),domClass1(1));
                                end
                            end 
                            if c1 < 10 %prev is full class - change 2nd class
                                if ismember(mixedClasses(k),mixedClassesTab(:,domClassFull(1))) == 0 
                                    if domClass2 == c1
                                        s2Probs(domClass2(1)) = 0;
                                        domClass2b = find(s2Probs == max(s2Probs));
                                        classImg(r,c,s2a:s2b) = findMixed(domClass2b(1),c1);
                                    else
                                        classImg(r,c,s2a:s2b) = findMixed(domClass2(1),c1);
                                    end
                                end
                            end
                        else %middle segment,
                            c1 = changeClasses(mixedChangeOrder(k)-1); 
                            c2 = changeClasses(mixedChangeOrder(k));
                            c3 = changeClasses(mixedChangeOrder(k)+1);
                            s1a = changeloc(mixedChangeOrder(k)-1);
                            s1b = changeloc(mixedChangeOrder(k))-1;
                            s2a = changeloc(mixedChangeOrder(k));
                            if length(changeloc) == mixedChangeOrder(k) %next segment is last segment
                                s2b = 39;
                            else
                                s2b = changeloc(mixedChangeOrder(k)+1)-1;
                            end
                            %if c1 and c3 are full class, then c2 should be
                            %mixed between them, but if full classes are
                            %the same, then just make sure it has perm
                            %class and 2nd most likely class for segment
                            %
                            if c1 < 10 &  c3 < 10 & c1 ~= c3
                                classImg(r,c,s2a:s2b) = findMixed(c1,c3);
                            elseif c1 < 10 & c1 == c3
                                s2Probs = mean((ccdcProbs(s2a:s2b,:)+anProbs(s2a:s2b,:))./2,1);
                                domClass2a = find(s2Probs == max(s2Probs));
                                s2Probs(domClass2a) = 0;
                                domClass2b = find(s2Probs == max(s2Probs));
                                if c1 ~= domClass2a(1)
                                     classImg(r,c,s2a:s2b) = findMixed(c1,domClass2a(1));
                                else
                                     classImg(r,c,s2a:s2b) = findMixed(c1,domClass2b(1));
                                end
                            elseif c1 < 10 & c3 > 10 %following class is mixed and will be corrected later, so just look at previous
                                s2Probs = mean((ccdcProbs(s2a:s2b,:)+anProbs(s2a:s2b,:))./2,1);
                                domClass2 = find(s2Probs == max(s2Probs));
                                if domClass2(1) == c1
                                    s2Probs(domClass2(1)) = 0;
                                    domClass2b = find(s2Probs == max(s2Probs));
                                    classImg(r,c,s2a:s2b) = findMixed(c1,domClass2b(1));
                                else
                                    classImg(r,c,s2a:s2b) = findMixed(c1,domClass2(1));
                                end
                            else %previous and following mixed, align with previous
                                s1Probs = mean((ccdcProbs(s1a:s1b,:)+anProbs(s1a:s1b,:))./2,1);
                                s2Probs = mean((ccdcProbs(s2a:s2b,:)+anProbs(s2a:s2b,:))./2,1);
                                segProbs = mean((ccdcProbs(s1a:s2b,:)+anProbs(s1a:s2b,:))./2,1);
                                domClassFull = find(segProbs == max(segProbs));
                                domClass1 = find(s1Probs == max(s1Probs));
                                domClass2 = find(s2Probs == max(s2Probs));
                                if ismember(mixedClasses(k),mixedClassesTab(:,domClassFull(1))) == 0 
                                    if domClassFull(1) == domClass2(1)
                                        s2Probs(domClass2)=0;
                                        domClass2b = find(s2Probs == max(s2Probs));
                                        classImg(r,c,s1a:s1b) = findMixed(domClass2b(1),domClassFull(1));
                                    else
                                        classImg(r,c,s1a:s1b) = findMixed(domClass2(1),domClassFull(1));
                                    end
                                end
                            end
                        end
                    end
                    %is it temp/perm?
                    if changeClasses(1) == changeClasses(end)
                        %gradual mixed change - temporary
                        changeType(r,c,1) = 3;
                    else
                        %gradual mixed change - permanent
                        changeType(r,c,1) = 4;
                    end
                end
            end
            %toc
        end
        %toc
    end
    toc
    save(['classImgSm-g',gridNums{gr},'.mat'],'classImg','changeType');
end
disp('done with smoothing')
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

%% - assigns mixed class based on 2 dominant classes
% uses position instead of class label
function mixClass = findMixed(c1,c2)
      if c1 == 1 & c2 == 2
          mixClass = 11; %water/emergent
      elseif c1 == 2 & c2 == 1
          mixClass = 11; %water/emergent
      elseif c1 == 1 & c2 == 3
          mixClass = 12; %water/forest-scrub
      elseif c1 == 3 & c2 == 1
          mixClass = 12; %water/forest-scrub
      elseif c1 == 1 & (c2 == 4 | c2 == 5)
          mixClass = 13; %water/other
      elseif c2 == 1 & (c1 == 4 | c1 == 5)
          mixClass = 13; %water/other
      elseif c1 == 2 & c2 == 3
          mixClass = 14;  %emergent/forest-scrub
      elseif c1 == 3 & c2 == 2
          mixClass = 14;  %emergent/forest-scrub
      elseif c1 == 2 & (c2 == 4 | c2 == 5)
          mixClass = 15; %emergent/other
      elseif c2 == 2 & (c1 == 4 | c1 == 5)
          mixClass = 15; %emergent/other
      elseif c1 == 3 & (c2 == 4 | c2 == 5)
          mixClass = 16; %forest-scrub/other
      elseif c2 == 3 & (c1 == 4 | c1 == 5)
          mixClass = 16; %forest-scrub/other
      elseif c1 == 4 & c2 == 5
          mixClass = 4; %not really mixed because 4 and 5 combined later
      elseif c1 == 5 & c2 == 4
          mixClass = 5; %not really mixed because 4 and 5 combined later
      end
end

%% function that gets ccdc and annual probabilities for a pixel


function [ccdcProbs, anProbs] = getProbs(ccdcPix, anPix)
    %build matrix
    ccdcProbs = zeros(39,5);
    anProbs = zeros(39,5);
    for yr = 1:39
        ccdcProbs(yr,:) = ccdcPix((yr-1)*6+2:yr*6); 
        anProbs(yr,:) = anPix((yr-1)*5+1:yr*5);
    end
    %fill missing values
    for k = 1:5
        ccdcProbs(:,k) = fillmissing(ccdcProbs(:,k),'previous');
        ccdcProbs(:,k) = fillmissing(ccdcProbs(:,k),'next');
        anProbs(:,k) = fillmissing(anProbs(:,k),'previous');
        anProbs(:,k) = fillmissing(anProbs(:,k),'next');
    end
end