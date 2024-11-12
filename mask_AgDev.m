%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to identify pixels that remain as class 4
% for the full time series - not interested in mapping these b/c never
% part of wetland area
% Courtney Di Vittorio
% Oct 15 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 
close all
%cd('E:\NASA_Coastal_Marshes\finalizeClassification_repSites\Barataria\Final_Maps\mosaicked_final')
cd('E:\Saeed\Coastwide_Classification\ws5_new\final_maps\mosaiced_final')
%get georeference so can save later for export
[tmpChange,R] = geotiffread('Mosaic_changeType.tif');
tmpImg = tmpChange;
%read in each year and stack
imgStack = zeros(size(tmpImg,1),size(tmpImg,2),39);
for yr = 1985:2023
    tmpImg = geotiffread(['Mosaic_',num2str(yr),'.tif']);
    imgStack(:,:,yr-1984) = tmpImg;
end

maskAgDev = zeros(size(tmpImg,1),size(tmpImg,2));
%make value 1 if sum of stack = 39*4, means that it was classified as other
%the entire 40 yrs, verified by a changeType = 0 to avoid possibility of
%other classes summing to same value
%maskAgDev(sum(imgStack,3) == (39*4)) = 1;
maskAgDev(sum(imgStack,3) == (39*4) & tmpChange == 0) = 1;
length(find(maskAgDev(:)==1))/length(maskAgDev(:)) %display percent of pixels masked

figure 
imagesc(maskAgDev)
colorbar

%%
cd('E:\Saeed\Coastwide_Classification\ws5_new\final_maps\')
%save as geotiff
geotiffwrite('maskAgDev.tif',maskAgDev,R)