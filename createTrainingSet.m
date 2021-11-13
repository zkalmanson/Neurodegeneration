function [idx,C] = createTrainingSet(k)
    %% Harvesting training set data  

    parent = pwd;
    image_names = dir(fullfile(parent,'*.tif'));
    image_names = natsortfiles({image_names.name});

    cropped = dendriteCropping(1); 
    [featureNum,~,~,~] = blebDetection(cropped);
    featuresMaster = zeros(10,featureNum);
    for i = 1:length(image_names)
        croppedImg = dendriteCropping(i); 
        [~,~,imgFeatures,~,~] = blebDetection(croppedImg); 
        featuresMaster = [featuresMaster;imgFeatures]; 

    end 
    featuresMaster = featuresMaster(all(featuresMaster,2),:);

    [idx,C] = kmeans(featuresMaster,k); 

end 
