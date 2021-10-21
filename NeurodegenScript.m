%% Get neurodegen data 

clear
clc
close all
parent = pwd;
image_names = dir(fullfile(parent,'*.tiff'));
image_names = natsortfiles({image_names.name});
blebCounts = zeros(length(image_names),1); 
dendHealth = zeros(length(image_names),1); 
saveLoc = 'C:\Users\zkalm\Desktop\DendriteDetection\Data\'; %CHANGE THIS TO WHATEVER FOLDER YOU WANT TO SAVE DATA IN! 
imgLoc = 'C:\Users\zkalm\Desktop\DendriteDetection\drive-download-20210723T234305Z-001'; %CHANGE THIS TO WHATEVER FOLDER YOU HAVE EVERYTHING IN!
addpath(imgLoc);

for i = 1:length(image_names)
    [blebCounts(i,1),dendHealth(i,1)] = ASC_NeuronCode_2021_v1_SurfKaze2(i); 
end 
    

%% Plot data 
close all

figure()
boxplot(blebCounts); 
ylabel('Number of blebs','fontsize',10)
title('Total Blebbing per Image','fontsize',10)
fig1 = get(groot,'CurrentFigure');
saveas(fig1, strcat(saveLoc, 'bleb count.tif')); 

figure()
boxplot(dendHealth); 
ylabel('Dendrite percent remaining','fontsize',10)
title('Dendrite Deterioration','fontsize',10)
fig2 = get(groot,'CurrentFigure');
saveas(fig2, strcat(saveLoc,'dendrite health.tif')); 