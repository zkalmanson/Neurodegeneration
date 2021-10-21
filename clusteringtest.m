    
    image = '2.tif'; 
    img = imread(image);
    [rows,columns] = size(img); 
    surfThresh = 1500;  
    surfThreshDisp = 50; 
    threshArray = zeros(surfThresh/surfThreshDisp,1); 
    blebCounts = zeros(surfThresh/surfThreshDisp,1); 


    surfPoints = detectSURFFeatures(img, 'MetricThreshold', surfThresh);
    [features, points] = extractFeatures(img, surfPoints); 
    surfXLoc1 = round(points.Location(:,1)); %Finds x,y coordinates of POI 
    surfYLoc1 = round(points.Location(:,2));  


    for i=1:surfThresh/surfThreshDisp
        currSurfThresh = surfThresh-(surfThreshDisp*i);
        surfPoints = detectSURFFeatures(img, 'MetricThreshold', currSurfThresh); 
        [sfeatures, surf_points] = extractFeatures(img, surfPoints); 

        surfXLoc = round(surf_points.Location(:,1)); %Finds x,y coordinates of POI 
        surfYLoc = round(surf_points.Location(:,2));  
        
        idx = kmeans(sfeatures,2); %Organize points into clusters 

        img1 = img(surfYLoc(idx==1),surfXLoc(idx==1)); 
        
        
        if mean(img(surfYLoc(idx==1),surfXLoc(idx==1)),'All') < mean(img(surfYLoc(idx==2),surfXLoc(idx==2)),'All')
            blebType = 2; 
        else 
            blebType = 1; 
        end
        
        
        blebTypes(i) = blebType;  
        totalXLoc = surfXLoc(idx==blebType);

        blebcount = length(totalXLoc); 
        blebCounts(i) = blebcount; %Final number of blebs 
        threshArray(i) = currSurfThresh;

        
        figure;
        imshow(img); 
        hold on; 
        plot(surfXLoc(idx==1),surfYLoc(idx==1),'r.','MarkerSize',8); 
        plot(surfXLoc(idx==2),surfYLoc(idx==2),'b.','MarkerSize',8); 
         
    end 

 
figure(); 
plot(threshArray,blebCounts);
xlabel('threshold');
ylabel('points'); 
 


 
