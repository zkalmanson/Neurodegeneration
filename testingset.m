%% Classify new images 
% This script takes trained clustering algorithm and sorts new points 
parent = pwd;
image_names = dir(fullfile(parent,'*.tif'));
image_names = natsortfiles({image_names.name});

cropped = cell(1,length(image_names));  %List of cropped dendrite 
featureMatrixes = cell(1,length(image_names)); %List of features used for sorting 
x = {}; % X and Y coords of feature points 
y = {}; 
for i = 1:length(image_names)
    cropped{i} = dendriteCropping(i); %Get cropped images 
    [~,~,featureMatrixes{i},x{i},y{i}] = blebDetection(cropped{i}); %Get points and features 

end 

k_value = 2; 
[idx,C] = createTrainingSet(k_value); %Train sorting algorithm 
idx_test = {}; 
for i = 1:length(image_names)
    [~,idx_test{i}] = pdist2(C,featureMatrixes{i},'euclidean','Smallest',1); %Sort test images 
end 


%% Visual example 

for i=1:length(image_names) 
    img_num = i; 
    figure;
    imshow(cropped{img_num}); 
    hold on; 
    plot(x{img_num}(idx_test{img_num}==1),y{img_num}(idx_test{img_num}==1),'r.','MarkerSize',15); 
    plot(x{img_num}(idx_test{img_num}==2),y{img_num}(idx_test{img_num}==2),'b.','MarkerSize',15); 
end 
