function [numOfFeatures,blebMask,features,final_x,final_y] = blebDetection(cropped_dendrite)
    % 1. Get SURF and KAZE points from cropped image 
    % 2. Intensity based gaussian filtering 
    % 3. Further filtering for points only in MSER regions of interest  
    % 4. Watershed segmentation and allocate one point per region 
    
    %% GET FEATURE DETECTOR POINTS 
    
    img = cropped_dendrite; 
    [rows,columns] = size(img);
    filt = fspecial('gaussian',2,0.25); %Gaussian filter 
    dendFilter = imfilter(img,filt);
    
    surfPoints = detectSURFFeatures(img, 'MetricThreshold', 250); %Finds points of interest 
    [~, surf_points] = extractFeatures(img, surfPoints); %Store information about POI
    
    kazePoints = detectKAZEFeatures(img, 'Threshold', 0.000001); %Finds points of interest 
    [~, kaze_points] = extractFeatures(img, kazePoints); %Store information about POI
    
    surfXLoc = round(surf_points.Location(:,1)); %Finds x,y coordinates of POI 
    surfYLoc = round(surf_points.Location(:,2));  
    kazeXLoc = round(kaze_points.Location(:,1)); %Finds x,y coordinates of POI 
    kazeYLoc = round(kaze_points.Location(:,2)); 
    
    totalXLoc = [surfXLoc;kazeXLoc]; %Combine surf and kaze points into one master list 
    totalYLoc = [surfYLoc;kazeYLoc]; 
    
    %{  
    figure(); 
    imshow(img);
    hold on; 
    plot(totalXLoc, totalYLoc,'x', 'LineWidth', .5, 'MarkerSize', 5, 'MarkerEdgeColor','r');
    %} 
     
    
    %% FILTERING THE POINTS 
    
    highcoords = zeros(100000,2);       
    highThresh = 0.75 * mean(dendFilter); %High pass filter to get rid of points we don't care about 
     
    for i = 1:length(totalXLoc) % Store high pass points 
        if dendFilter(totalYLoc(i),totalXLoc(i)) > highThresh  % Store high pass points  
           highcoords(i,1) = totalXLoc(i); 
           highcoords(i,2) = totalYLoc(i); 
        end 
    end
    
    highcoords = highcoords(all(highcoords,2),:); %Removes zeros 
    
    
    %{
    figure();
    imshow(img);
    hold on; 
    plot(highcoords(:,1), highcoords(:,2),'x', 'LineWidth', .5, 'MarkerSize', 5, 'MarkerEdgeColor','y');
    %}   
    
    [regions,mserCC] = detectMSERFeatures(img,'ThresholdDelta',0.25); % Construct MSER regions of interest 
    
    stats = regionprops('table',mserCC,'All'); 
    circIdx = find([stats.Circularity] > 0.2); %Only want circular regions 
    circregions = regions(circIdx); 
    circpix = circregions.PixelList; 
    circpixlist = cell2mat(circpix); 
    circpixlistX = circpixlist(:,1); 
    circpixlistY = circpixlist(:,2); 

    MSERmask = zeros(rows,columns); 
    for i = 1:length(circpixlistX)
        MSERmask(circpixlistY(i),circpixlistX(i)) = 1; %Creating mser mask 
    end 
    MSERmask = logical(MSERmask); 
    MSERmaskCC = bwconncomp(MSERmask); 
    
    
    MSERlabel = labelmatrix(MSERmaskCC);
    MSERcoords = zeros(100000,2);  
    for i=1:length(highcoords(:,1))
        if MSERlabel(highcoords(i,2),highcoords(i,1)) ~= 0 % Remove points not in a MSER region mask
            MSERcoords(i,1) = highcoords(i,1); 
            MSERcoords(i,2) = highcoords(i,2); 
        end
    end 
    MSERcoords = MSERcoords(all(MSERcoords,2),:);

    
    filteredPoints(:,1) = MSERcoords(:,1);
    filteredPoints(:,2) = MSERcoords(:,2);
    
    %{
    figure();
    imshow(img);
    hold on; 
    plot(circregions,'showPixelList',true,'showEllipses',false)
    plot(filteredPoints(:,1), filteredPoints(:,2),'x', 'LineWidth', .5, 'MarkerSize', 5, 'MarkerEdgeColor','r');
    %} 
    
    
    %% Watershed segmention to allocate one point per "region" of dendrite 
    
    t = adaptthresh(img,0.15);
    dendBW = imbinarize(img,t);
    flatBW = bwareaopen(dendBW,10); 

    D = bwdist(~flatBW);
    D = -D;

    minMask = imextendedmin(D,0.85);
    D2 = imimposemin(D,minMask); 
    lblBW = watershed(D2); 
    lblBW(~flatBW) = 0; 
    labeled = bwlabel(lblBW); %Labels watershed regions
    
    for i=1:length(filteredPoints(:,1))
        if labeled(filteredPoints(i,2),filteredPoints(i,1)) ~= 0
            filteredPoints(i,4) = labeled(filteredPoints(i,2),filteredPoints(i,1)); %Displays what watershed region each coordinate is in 
        end
    end
    
    % WANT TO IMPROVE THIS SO IT TAKES "STRONGEST" POINT IN EACH REGION  
    filteredIndices = zeros(1000,1);
    for i=1:max(labeled,[],'all') % For each watershed region with multiple points, deletes all but first one. There is probably a better way to do this... 
        if length(find(filteredPoints(:,4) == i)) > 0
            tempPoints = find(filteredPoints(:,4) == i); %Temporary list to store indices of points in this region 
            filteredIndices(i,1) = tempPoints(1,1); %Take the first point and store index into new list 
        else
            filteredIndices(i,1) = 0; 
        end 
    end 
    
    filteredIndices = filteredIndices(all(filteredIndices,2),:); % Linear indices of final points 
    [row,col] = ind2sub([544 4],filteredIndices); % Convert back into matrix coords 
    
    final_x = filteredPoints(row,1); % Final coords 
    final_y = filteredPoints(row,2); 
    
    
    figure();
    imshow(img);
    hold on; 
    plot(final_x,final_y,'x', 'LineWidth', .5, 'MarkerSize', 5, 'MarkerEdgeColor','r');
    
    
    rgb = label2rgb(lblBW,'jet','w','shuffle');
    figure();
    alpha = imshow(rgb); 
    alpha.AlphaData = 0.3;
    hold on; 
    plot(final_x, final_y,'x', 'LineWidth', .5, 'MarkerSize', 5, 'MarkerEdgeColor','b');
    
    
    
    
    %% CREATE FINAL MASK WITH DETECTED POINTS
    
    blebCount = length(final_x); %Final number of points  
    
    %{
    figure();   
    imshow(img);
    hold on; 
    plot(final_x+5, final_y,'_', 'LineWidth', .5, 'MarkerSize', 3, 'MarkerEdgeColor','g');
    %} 
    
    
    blebMask = zeros(rows,columns); 
    for i = 1:length(final_x)
        blebMask(final_y(i),final_x(i)) = 1; 
    end 
    
    blebMask = logical(blebMask); % Final mask where all blebs are denoted 1, other pixels given a -11   
    
    
    %% K CLUSTERING TRAINING DATA, THIS SECTION IS WORK IN PROGRESS!! 
    
    for i=1:length(final_x)
        pixels(i) = img(final_y(i),final_x(i)); 
    end 
    
    numOfFeatures = 13; 
    features = zeros(length(pixels),numOfFeatures,'double'); 
    
    features(:,1) = pixels; 
    
    f1 = fspecial('gaussian',15,3); 
    f1img = imfilter(img,f1); 
    for i=1:length(final_x)
        features(i,2) = f1img(final_y(i),final_x(i)); 
    end 
    
    f2 = fspecial('gaussian',15,6); 
    f2img = imfilter(img,f2); 
    for i=1:length(final_x)
        features(i,3) = f2img(final_y(i),final_x(i)); 
    end 
    
    f3 = fspecial('gaussian',15,10); 
    f3img = imfilter(img,f3); 
    for i=1:length(final_x)
        features(i,4) = f3img(final_y(i),final_x(i)); 
    end 
    
    f4img = stdfilt(img); 
    for i=1:length(final_x)
        features(i,5) = f4img(final_y(i),final_x(i)); 
    end 
    
    f5 = fspecial('disk',5);  
    f5img = imfilter(img,f5);
    for i=1:length(final_x)
        features(i,6) = f5img(final_y(i),final_x(i)); 
    end 
 
    f6 = fspecial('sobel'); 
    f6img = imfilter(img,f6);
    for i=1:length(final_x)
        features(i,7) = f6img(final_y(i),final_x(i)); 
    end 
    
    f7 = entropyfilt(img);
    for i=1:length(final_x)
        features(i,8) = f7(final_y(i),final_x(i)); 
    end 
    
    
    f8 = imgaborfilt(img,2,0); %Features 9-14 are gabor filters with varying wavelength and phase 
    for i=1:length(final_x)
        features(i,9) = f8(final_y(i),final_x(i)); 
    end 
    
    f9 = imgaborfilt(img,2,90); 
    for i=1:length(final_x)
        features(i,10) = f9(final_y(i),final_x(i)); 
    end 
    
    f10 = imgaborfilt(img,2,45); 
    for i=1:length(final_x)
        features(i,11) = f10(final_y(i),final_x(i)); 
    end 
    
    f11 = imgaborfilt(img,2,180); 
    for i=1:length(final_x)
        features(i,12) = f11(final_y(i),final_x(i)); 
    end 
    
    f12 = imgaborfilt(img,2,120); 
    for i=1:length(final_x)
        features(i,13) = f12(final_y(i),final_x(i)); 
    end 
    
    %{
    idx = kmeans(features,3);
    
    figure;
    imshow(img); 
    hold on; 
    plot(final_x(idx==1),final_y(idx==1),'r.','MarkerSize',15); 
    plot(final_x(idx==2),final_y(idx==2),'b.','MarkerSize',15); 
    plot(final_x(idx==3),final_y(idx==3),'c.','MarkerSize',15); 
    %} 
    
end

    
