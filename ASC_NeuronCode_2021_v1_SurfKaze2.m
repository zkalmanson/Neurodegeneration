% Attempt code on new images
% September 16, 2020
% v11 - attempt to use neighborhood binning

function [finalBlebCount, brk] = ASC_NeuronCode_2021_v1_SurfKaze2(img_number)
    %% GET IMAGES
    close all
    %% Load images and make grayscale for processing
    parent = pwd;
    image_names = dir(fullfile(parent,'*.tiff'));
    image_names = natsortfiles({image_names.name});
    img_num = img_number; %which image to look at
    im_ref_num = 7; %match hist eq to img 1
    split = 20;
    meansplit = 20;
    meanrange = 40;
    changeval = 10;
    minprom = 15;
    file=fullfile(parent,image_names{img_num});
    im = (imread(file));
    img = im(:,:,2);
    %imshow(img)
    
    %% Binarize image and rotate to be vertical
    close all
    imbw = imbinarize(img,0.1); %make high to only get brights section 
    close all
    prps = regionprops(imbw,'Area','Orientation'); % Obtain orientation so can rotate to make vertical
    img_prps(:,1) = [prps.Area];
    img_prps(:,2) = [prps.Orientation];
    [max_area, max_idx] = max(img_prps(:,1));
    max_agl = img_prps(max_idx,2);
    imrot  = imrotate(img,90-max_agl,'nearest');
    %imshow(imrot);
    
    %% Cut image into 2 - get side with dendrites
    close all
    % Smaller binarization for center of cell body
    imbw2 = imbinarize(imrot,0.9); 

    % Larger binarization to get larger Bounding Box (used to crop images)
    close all
    imbw3 = imclearborder(imbinarize(imrot,0.4));
    % imshow(imbw3)
    
    % Get regionprops data for both binary images
    rot_prps = regionprops(bwareaopen(imbw2,10),'all');
    rot_area(:,1) = [rot_prps.Area]; % puts into array
    [max_a2, max_idx2] = max(rot_area(:,1)); % Finds max object (will be the cell body)
    rot_cell = struct2cell(rot_prps); % makes into array
    rot_bb = cell2mat(rot_cell(3,max_idx2)); % Finds the width of the bounding box of the largest object
    rot_center = cell2mat(rot_cell(2,max_idx2)); % Finds the center y value of the bounding box
    
    % Finds larger binarized object bounding box for cropping purposes
    rot_bb_prps = regionprops(bwareaopen(imbw3,10),'Area','BoundingBox');
    rot3_cell = struct2cell(rot_bb_prps);
    rot2_area(:,1) = rot3_cell(1,:);
    [max_a3, max_idx3] = max(cell2mat(rot2_area(:,1)));
    x_rot_bb = cell2mat(rot3_cell(2,max_idx3));

    % Cropping locations from regionprops
    center_y = rot_center(2);
    x_range = [rot_center(1)-1.5*x_rot_bb(3) rot_center(1)+1.5*x_rot_bb(3)]; %uses 1.5x the x length of the boudning box
    y_range = [rot_bb(2) rot_bb(2)+1*rot_bb(4)]; %determines the bounding box of the cell body
    x_mag = x_range(2)-x_range(1);
    top_y = 1; %first 'row'
    bot_y = length(imrot(:,1)); % Last 'row' of the image

    % Now create 2 images
    imtop = imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]);
    imbot = imcrop(imrot, [x_range(1) center_y x_mag bot_y-center_y]);  %use y range to get rid of cell body/ center_y to include cell body
    %imshowpair(imtop,imbot,'montage')

    %% Decide which side has dendrites
    % Binarize both top and bottom images
    close all
    T_top = adaptthresh(imtop(5:size(imtop,1)-4,5:size(imtop,2)-4),0.22,'neighborhoodsize',[31 31]);
    T_bot = adaptthresh(imbot(5:size(imbot,1)-4,5:size(imbot,2)-4),0.22,'neighborhoodsize',[31 31]);
    top_bw = imclearborder(bwareaopen(imbinarize(imtop(5:size(imtop,1)-4,5:size(imtop,2)-4),T_top),250));
    bot_bw = imclearborder(bwareaopen(imbinarize(imbot(5:size(imbot,1)-4,5:size(imbot,2)-4),T_bot),250));
    %imshowpair(top_bw,bot_bw,'montage')
    
    % Use regionprops data to determine which side is correct
    close all
    top_prps = regionprops(top_bw,'all');
    bot_prps = regionprops(bot_bw,'all');
    
    if isempty(top_prps) % Nothing is binarized
        top_stats(:,1) = 0;
        top_stats(:,2) = 0;
    else % Get area and orientation
        top_stats(:,1) = [top_prps.Area];
        top_stats(:,2) = [top_prps.Orientation];
    end
   
    if isempty(bot_prps) % Nothing binarized
        bot_stats(:,1) = 0;
        bot_stats(:,2) = 0;
    else % Get area and orientation
        bot_stats(:,1) = [bot_prps.Area];
        bot_stats(:,2) = [bot_prps.Orientation];
    end

    % Finds largest 4 objects (or less)
    [top_a_max, top_idx] = maxk(top_stats(:,1),4);
    [bot_a_max, bot_idx] = maxk(bot_stats(:,1),4);
    
    % Gets orientation of largest 4 objects
    top_ang = top_stats(top_idx,2);
    bot_ang = bot_stats(bot_idx,2);
    
    % Obtains standard deviations
    top_std = std(abs(top_ang));
    bot_std = std(abs(bot_ang));
    
    % Size of each images row x col
    topSize = size(top_bw,1);
    botSize = size(bot_bw,1);
    
    % Proportional size to original image
    topSizeProp = topSize/(topSize+botSize);
    botSizeProp = botSize/(topSize+botSize);
    
    % Complete if statements to determine which is correct size
    % If nothing binarized, other side must be the dendrites
    if size(top_prps,1) == 0 && mean(abs(bot_ang)) >= 75
        imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
        im_dend = imbot;
        bordersToKeep = [1 0 0 0];
    elseif size(bot_prps,1) == 0 && mean(abs(top_ang)) >= 75
        imtop = imrotate(imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
        im_dend = imtop;
        bordersToKeep = [1 0 0 0];
    elseif botSizeProp > 0.6 && ((size(top_prps,1) == 0 || size(bot_prps,1) == 0))
        imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
        im_dend = imbot;
        bordersToKeep = [1 0 0 0];
    elseif topSizeProp > 0.6 && ((size(top_prps,1) == 0 || size(bot_prps,1) == 0))
        imtop = imrotate(imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
        im_dend = imtop;
        bordersToKeep = [1 0 0 0];
    elseif size(top_prps,1) == 1 && size(top_prps,1)<size(bot_prps,1) && mean(abs(bot_ang)) >= 75
        imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
        im_dend = imbot;
        bordersToKeep = [1 0 0 0];        
    elseif size(bot_prps,1) == 1 && size(top_prps,1)>size(bot_prps,1) && mean(abs(top_ang)) >= 75
        imtop = imrotate(imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
        im_dend = imtop;
        bordersToKeep = [1 0 0 0];   
    else
        % Determine which image to use by which has the most large objects closest
        % to parallel
        bordersToKeep = [1 0 0 0];
        if top_std < bot_std && mean(abs(top_ang)) >= 80
            imtop = imrotate(imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
            im_dend = imtop;
        elseif top_std > bot_std && mean(abs(bot_ang)) >= 80
            imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
            im_dend = imbot;
        elseif abs(mean(abs(top_ang))-90)  < abs(mean(abs(bot_ang)) - 90)
            imtop = imrotate(imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
            im_dend = imtop;
        elseif abs(mean(abs(top_ang))-90)  > abs(mean(abs(bot_ang)) - 90)
            imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
            im_dend = imbot;
        else
            %imshowpair(top_bw,bot_bw,'montage');
            imgprompt = ('Can not pick out which side has dendrites! \nPlease pick which side has dendrites: Left (1) or Right (2)');
            user_img = input(imgprompt);
            user_img = str2double(user_img);
            while user_img ~= 1 || user_img ~= 2
                if user_img ~= 1 || user_img ~= 2
                    break
                else
                    user_img = input(imgprompt);
                    user_img = str2double(user_img); 
                end
            end
            if user_img == 1
                imtop = imrotate(imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
                im_dend = imtop; 
            elseif user_img == 2
                imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
                im_dend = imbot;            
            end
        end
    end
    close all
    %imshow(im_dend)
    
    %{ 
    % User checks if it chose the right size
    userprompt = ('Is this the correct picture? [Y]/N \n');
    usercheck = input(userprompt,'s');
    usercheck = lower(usercheck);
    
    % Do not need to enter 'y' - can just hit enter
    if isempty(usercheck)
        usercheck = 'y';
    end
    while usercheck ~= 'y' || usercheck ~= 'n'
        if usercheck ~= 'y' || usercheck ~= 'n'
            break
        else
        userprompt = ('Is this the correct picture? (Y/N) \n');
        usercheck = input(userprompt);
        usercheck = lower(usercheck);
        end
    end
    topcheck = mean2(imtop);
    botcheck = mean2(imbot);
    dendcheck = mean2(im_dend);
    if usercheck == 'y'
        close all
    elseif usercheck == 'n' && sum(topcheck-dendcheck) == 0
        close all
        imbot = imcrop(imrot, [x_range(1) y_range(2) x_mag bot_y-y_range(2)]);
        im_dend = imbot;       
    elseif usercheck == 'n' && sum(botcheck-dendcheck) == 0
        close all
        imtop = imrotate(imcrop(imrot,[x_range(1) top_y x_mag y_range(1)-top_y]),180);
        im_dend = imtop;

    end
    %} 
    %% Now we have the right part of the image - time to crop down more
    % Get rid of cell body before next round of cropping 
    close all
    
    % Binarize image to obtain only remainder of cell body (2 different ways)
    imbw4 = imbinarize(im_dend,0.2);
    imbw5 = imerode(imbw4,strel('rectangle',[2,10]));
    imbw5clear = imclearborder(imbw5); % Clears border
    
    % Remove only top border
    imbw4edgepad = padarray(imbw4,[bordersToKeep(1), bordersToKeep(4)],'pre');
    imbw4edgepad = padarray(imbw4edgepad,[bordersToKeep(3), bordersToKeep(2)],'post');
    imbw4edgepad2 = imclearborder(imbw4edgepad);
    topLeft = [bordersToKeep(1), bordersToKeep(4)]+1;
    botRight = topLeft + [size(imbw4,1), size(imbw4,2)] - 1;
    imbw4clear =imbw4edgepad2(topLeft(1):botRight(1), topLeft(2):botRight(2)); 
    imbw4clear = bwareaopen(imbw4clear,5); % Removes small extra binarizations
    %imshow(imbw4clear)
    
    % Used to remove cell body
    imbw5edge = imbw5;
    imbw5edge(imbw5clear) = 0; % Pixels that do not touch the border of imbw5 become zero
    imbw5edge = imerode(imbw5edge,strel('line',20,0));
    %imshowpair(imbw4,imbw5edge,'montage')

    % Clears border except top of image
    imedgepad = padarray(imbw5edge,[bordersToKeep(1), bordersToKeep(4)],'pre');
    imedgepad = padarray(imedgepad,[bordersToKeep(3), bordersToKeep(2)],'post');
    imedgepad2 = imclearborder(imedgepad);
    %imshowpair(imbw5edge,imedgepad,'montage')
    topLeft = [bordersToKeep(1), bordersToKeep(4)]+1;
    botRight = topLeft + [size(imbw5edge,1), size(imbw5edge,2)] - 1;
    imbw5edge =imedgepad2(topLeft(1):botRight(1), topLeft(2):botRight(2));
    %imshow(imbw5edge)
    
    % Finds the x boundaries of where it needs to be cropped 
    clearmax = zeros(size(imbw4clear,2),1); % Preallocates to find x indicies 
    for ii = 1:size(imbw4clear,2)
       clearmax(ii) = max(imbw4clear(:,ii)); % Uses original imbw4
    end
    clearidx(1) = find(clearmax,1);
    clearidx(2) = find(clearmax,1,'last'); % will be used later
    
    % Obtain bounding box to of cell body to be used for cropping
    close all
    imbw5props = regionprops(imbw5edge,'Area','BoundingBox');
    rot5_cell = struct2cell(imbw5props);
    rot5_area(:,1) = rot5_cell(1,:);
    [max_a5, max_idx5] = max(cell2mat(rot5_area(:,1)));
    cellbody_bb = cell2mat(rot5_cell(2,max_idx5));
    imbw6 = imcrop(imbw4,[clearidx(1)-10, cellbody_bb(2)+cellbody_bb(4)-10, clearidx(2)-clearidx(1)+20, size(imbw4,1)]); 
    % Only dendrites should be in imbw6
    %imshowpair(imbw4,imbw6,'montage')

    % Dilate dendrites and remove objects that touch any border besides the top
    im_dilpad = padarray(imbw6,[bordersToKeep(1), bordersToKeep(4)],'pre');
    im_dilpad = padarray(im_dilpad,[bordersToKeep(3), bordersToKeep(2)],'post');
    im_dil2 = imclearborder(im_dilpad);
    %imshowpair(imbw6,im_dilpad,'montage')
    topLeft = [bordersToKeep(1), bordersToKeep(4)]+1;
    botRight = topLeft + [size(imbw6,1), size(imbw6,2)] - 1;
    im_dil =im_dil2(topLeft(1):botRight(1), topLeft(2):botRight(2));
    %imshow(im_dil)
    
   
    
    %imshowpair(imbw5,imbw4clear,'montage')    %im_dend = imcrop(im_dend,[cellbody_bb(1)-30, cellbody_bb(2)+cellbody_bb(4)-10, cellbody_bb(3)+60, size(imbw4,1)]); 
    
    
    % Crops original image to focus on denbdrites - again
    im_dend = imcrop(im_dend,[clearidx(1)-10, cellbody_bb(2)+cellbody_bb(4)-10, clearidx(2)-clearidx(1)+20, size(imbw4,1)]); 
    
    % Closes image to get bounding box + orientation for further cropping
    im_close = imclose(im_dil,strel('disk',100,8));
    %figure()
    %imshowpair(im_dend,im_close)

    % Get region props data and rotate iomage again to get vertical
    close all
    close_or = regionprops(bwareaopen(im_close,5),'orientation');
    close_area = regionprops(bwareaopen(im_close,5),'area');
    [close_a, close_idx] = max(cell2mat(struct2cell(close_area)));
    close_o = cell2mat(struct2cell(close_or(close_idx)));
    close_sign = sign(close_o);
    if abs(90-abs(close_o)) > 5
        im_close = imrotate(im_close,close_sign*(90-abs(close_o)));
        im_dend = imrotate(im_dend,close_sign*(90-abs(close_o)));
        %imshowpair(im_dend,im_close);
    end
    % Bounding box collection
    close_props = regionprops(bwareaopen(im_close,5),'boundingbox');
    close_bb = cell2mat(struct2cell(close_props(close_idx)));
    close_bb(1) = close_bb(1)-10; % add a little padding in x
    close_bb(3) = close_bb(3)+20; 
    close_bb(4) = close_bb(4)-10;

    % Crop to final image for processing
    dend_crop = imcrop(im_dend,close_bb);
    %imshow(dend_crop)
    
    %% Implement Bleb Detection 
   
    [rows,columns] = size(dend_crop);
    filt = fspecial('gaussian',5,0.75); %Filter for finding blebs
    dendFilt = imfilter(dend_crop,filt);
    
    surfPoints = detectSURFFeatures(dend_crop, 'MetricThreshold', 1000); %Finds points of interest 
    [~, surf_points] = extractFeatures(dend_crop, surfPoints); %Store information about POI
    
    kazePoints = detectKAZEFeatures(dend_crop); %Finds points of interest 
    [~, kaze_points] = extractFeatures(dend_crop, kazePoints); %Store information about POI
    
    surfXLoc = round(surf_points.Location(:,1)); %Finds x,y coordinates of POI 
    surfYLoc = round(surf_points.Location(:,2));  
    kazeXLoc = round(kaze_points.Location(:,1)); %Finds x,y coordinates of POI 
    kazeYLoc = round(kaze_points.Location(:,2)); 
    
    totalXLoc = [surfXLoc;kazeXLoc]; 
    totalYLoc = [surfYLoc;kazeYLoc]; 
    
    %{  
    figure();
    imshow(img);
    hold on; 
    plot(totalXLoc, totalYLoc,'x', 'LineWidth', .5, 'MarkerSize', 5, 'MarkerEdgeColor','r');
    %} 
   
    
    lowcoords = zeros(100000,2); 
    lowLength = 0; 
     
    highcoords = zeros(100000,2);  
    highLength = 0; 
     
    highThresh = 2 * mean(dendFilt); 
     
    for i = 1:length(totalXLoc) % Store low pass points 
        if dendFilt(totalYLoc(i),totalXLoc(i)) < 150 && dendFilt(totalYLoc(i),totalXLoc(i)) > 100 %120   
           lowcoords(i,1) = totalXLoc(i); 
           lowcoords(i,2) = totalYLoc(i); 
           lowLength = lowLength + 1; 
        end 
        if dendFilt(totalYLoc(i),totalXLoc(i)) > highThresh  % Store high pass points (170) (200) 
           highcoords(i,1) = totalXLoc(i); 
           highcoords(i,2) = totalYLoc(i); 
           highLength = highLength + 1; 
        end 
    end
     
    lowcoords = lowcoords(all(lowcoords,2),:);
    highcoords = highcoords(all(highcoords,2),:);
    
    filter_X = [lowcoords(:,1);highcoords(:,1)]; 
    filter_y = [lowcoords(:,2);highcoords(:,2)]; 
    
    %{
    figure();
    imshow(img);
    hold on; 
    plot(highcoords(:,1), highcoords(:,2),'x', 'LineWidth', .5, 'MarkerSize', 5, 'MarkerEdgeColor','y');
    %} 
    
   
    
    [regions,mserCC] = detectMSERFeatures(dend_crop,'ThresholdDelta',1); % Construct MSER regions 
    
    stats = regionprops('table',mserCC,'All'); 
    circIdx = find([stats.Circularity] > 0.3); 
    circregions = regions(circIdx); 
    circpix = circregions.PixelList; 
    pixlist = cell2mat(circpix); 
    pixlistX = pixlist(:,1); 
    pixlistY = pixlist(:,2); 

    
    mask = zeros(rows,columns); 
    for i = 1:length(pixlistX)
        mask(pixlistY(i),pixlistX(i)) = 1; 
    end 
    mask = logical(mask); 
    maskCC = bwconncomp(mask); 
    
    
    mserlabel = labelmatrix(maskCC);
    msercoords = zeros(100000,2);  
    for i=1:length(highcoords(:,1))
        if mserlabel(highcoords(i,2),highcoords(i,1)) ~= 0 % Remove points not in a MSER region 
            msercoords(i,1) = highcoords(i,1); 
            msercoords(i,2) = highcoords(i,2); 
        end
    end 
    msercoords = msercoords(all(msercoords,2),:);

    
    lowHigh(:,1) = msercoords(:,1);
    lowHigh(:,2) = msercoords(:,2);
    
    %{ 
    figure();
    imshow(img);
    hold on; 
    plot(circregions,'showPixelList',true,'showEllipses',false)
    plot(lowHigh(:,1), lowHigh(:,2),'x', 'LineWidth', .5, 'MarkerSize', 5, 'MarkerEdgeColor','r');
    %} 
    
    
    t = adaptthresh(dend_crop,0.5);
    dendBW = imbinarize(dend_crop,t);
    flatBW = bwareaopen(dendBW,10); 

    D = bwdist(~flatBW);
    D = -D;

    minMask = imextendedmin(D,0.85);
    D2 = imimposemin(D,minMask); 
    lblBW = watershed(D2); 
    lblBW(~flatBW) = 0; 
    labeled = bwlabel(lblBW); %Labels watershed regions
    
    for i=1:length(lowHigh(:,1))
        if labeled(lowHigh(i,2),lowHigh(i,1)) ~= 0
            lowHigh(i,4) = labeled(lowHigh(i,2),lowHigh(i,1)); %Displays what watershed region each coordinate is in 
        end
    end 

    [~, w] = unique(lowHigh(:,4), 'stable' );
    duplicate_indices = setdiff(1:numel(lowHigh(:,4)), w ); %Creates array of duplicate coordinate indices 

    lowHigh(duplicate_indices,:) = []; %Removes duplicate coordinates 
    
    %{ 
    rgb = label2rgb(lblBW,'jet','w','shuffle');
    figure();
    alpha = imshow(rgb); 
    alpha.AlphaData = 0.3;
    hold on; 
    plot(lowHigh(:,1), lowHigh(:,2),'x', 'LineWidth', .5, 'MarkerSize', 5, 'MarkerEdgeColor','b');
    %} 
    
    
    
    x_locFinal = lowHigh(:,1); 
    y_locFinal = lowHigh(:,2); 
    
    blebCount = length(x_locFinal); %Final number of blebs 
    
    %{
    figure();   
    imshow(dend_crop);
    hold on; 
    plot(x_locFinal+5, y_locFinal,'_', 'LineWidth', .5, 'MarkerSize', 3, 'MarkerEdgeColor','g');
    %} 
    
    blebMask = zeros(rows,columns); 
    for i = 1:length(x_locFinal)
        blebMask(y_locFinal(i),x_locFinal(i)) = 1; 
    end 
    
    blebMask = logical(blebMask); % Final mask where all blebs are denoted 1, other pixels given a -11    
    

    %% Find 4 max points for each 'y-value' (row)
    
    close all
    % Preallocation
    max_pt = zeros(size(dend_crop));
    max_loc = zeros(size(dend_crop,1),4);
    max_pic = zeros(size(dend_crop));
    
    % Find max intensity value at each row
    for ii = 1:size(dend_crop,1)
        max_pt(ii,:) = islocalmax(dend_crop(ii,:),'maxnumextrema',4,'minseparation',5, ...
            'minprominence',minprom); % Max 4 points
        length_max = length(find(max_pt(ii,:)== 1)); % counts number of points
        length_max_up = 5-length_max;

        if length_max == 0 % if nothing found - put max points at first row
            max_loc(ii,:) = [1 1 1 1];
            
        elseif length_max < 4 % if under 4 found - put missing at row 1
            %gets locations of possible dendrites
            max_loc(ii,1:4-length_max) = 1;
            max_loc(ii,length_max_up:4) = find(max_pt(ii,:) == 1);

        else % if 4 found - index for binary image is these 4
            max_loc(ii,:) = find(max_pt(ii,:) == 1);
        end
        
        % Makes max index points = 1 and other = 0
        max_pic(ii,max_loc(ii,:)) = 1;
    end
    
    % Create binary image and removes 'objects' under 2 pixels
    max_pic = bwareaopen(logical(max_pic(:,2:size(max_pic,2))),2); 
    %imshowpair(dend_crop,max_pic)

    %% Attempt at placing index points in correct column for individual dendrite
    % Creates new index variables to not mess with original
    max_loc2 = max_loc;
    max_loc3 = max_loc;

    max_loc2(any(max_loc2==1,2),:) = []; % Removes points where max_idx = 1
    max_loc_avg = round(median(max_loc2,1)); % Gets the median location for each of the 4 cols
    
    % Obtains average locations in between dendrites
    dend_loc_avg(1) = round(mean(max_loc_avg(1:2))); 
    dend_loc_avg(2) = round(mean(max_loc_avg(2:3)));
    dend_loc_avg(3) = round(mean(max_loc_avg(3:4)));
    
    % Creates column ranges to determine where dendrites should be placed
    max_avg_rng(1,:) = [2 dend_loc_avg(1)];
    max_avg_rng(2,:) = [dend_loc_avg(1) dend_loc_avg(2)];
    max_avg_rng(3,:) = [dend_loc_avg(2) dend_loc_avg(3)];
    max_avg_rng(4,:) = [dend_loc_avg(3) size(max_pic,2)];

    max_pic2 = zeros(size(dend_crop)); % Preallocates to make image
    q = 0;
    A = 0;
    Q = size(max_loc2,1);
    obsRng = 7;
    nuer_range = zeros(size(max_pic,1),2,4);
    nuer_idx = zeros(size(max_loc));
    
    % For loop to place pixels in correct range
    for kk = 1:size(max_pic,1)

        if range(max_loc3(kk,:)) == 0 % If all idx = 1 (no pt max pt detected), use average pont for this y value
            for mm = 1:4
                nuer_range(kk,:,mm) = [max_loc_avg(mm)-obsRng max_loc_avg(mm)+obsRng];
            end

        elseif min(max_loc3(kk,:)) == 1 % if less than 4 points detected

            % attempts to place max location in correct vector index when less than
            % 4 nuerons detected      

            for nn = 1:4
               loc_max_idx = max_loc(kk,nn); % Creates temp variable for this row
               
               for oo = 1:4
                   if loc_max_idx < max_avg_rng(oo,2) && loc_max_idx >= max_avg_rng(oo,1)

                       if nuer_idx(kk,oo) == 0
                            nuer_idx(kk,oo) = loc_max_idx;
                       elseif nuer_idx(kk,oo) > loc_max_idx
                           nuer_idx(kk,oo-1) = loc_max_idx;
                       elseif nuer_idx(kk,oo) < loc_max_idx && oo == 1 && nuer_idx(kk,oo+1) == 0
                           nuer_idx(kk,oo+1) = loc_max_idx;
                       elseif nuer_idx(kk,oo) < loc_max_idx && oo == 1 && nuer_idx(kk,oo+1) ~= 0
                           if nuer_idx(kk,oo+2) == 0
                               nuer_idx(kk,oo+2) = loc_max_idx;
                           elseif nuer_idx(kk,oo+3) == 0
                               nuer_idx(kk,oo+3) = loc_max_idx;   
                           end
                       elseif nuer_idx(kk,oo) < loc_max_idx && nuer_idx(kk,oo-1) == 0
                           nuer_idx(kk,oo-1) = nuer_idx(kk,oo);
                           nuer_idx(kk,oo) = loc_max_idx;
                       elseif nuer_idx(kk,oo) < loc_max_idx && nuer_idx(kk,oo-1) ~= 0 && oo > 2
                           nuer_idx(kk,oo-2) = nuer_idx(kk,oo-1);
                           nuer_idx(kk,oo-1) = nuer_idx(kk,oo);
                           nuer_idx(kk,oo) = loc_max_idx;
                       elseif nuer_idx(kk,oo) < loc_max_idx && nuer_idx(kk,oo-1) ~= 0 && oo <= 2 && nuer_idx(kk,oo+1) == 0
                           nuer_idx(kk,oo+1) = loc_max_idx;
                       elseif nuer_idx(kk,oo) < loc_max_idx && nuer_idx(kk,oo-1) ~= 0 && oo <= 2 && nuer_idx(kk,oo+1) ~= 0
                           if nuer_idx(kk,oo+2) == 0
                           nuer_idx(kk,oo+2) = loc_max_idx;
                           elseif nuer_idx(kk,oo+2) ~= 0 && oo == 1
                               nuer_idx(kk,oo+3) = loc_max_idx;
                           end
                       else
                           fprintf('ERROR: Please check max location rearrangement! \n');
                       end
                   end
               end
            end
        else
            nuer_idx(kk,:) = max_loc(kk,:);
        end
    end
    
    % Get average of neur_idx over non zero points
    neur_idx_avg = round([sum(nuer_idx(:,1))./nnz(nuer_idx(:,1)), sum(nuer_idx(:,2))./nnz(nuer_idx(:,2)), ...
        sum(nuer_idx(:,3))./nnz(nuer_idx(:,3)), sum(nuer_idx(:,4))./nnz(nuer_idx(:,4))]);

    %Attempt to show nuer_idx location with image
    nuer_idx2 = nuer_idx;
    nuer_idx2(nuer_idx2<=0) = 1;
    nuer_idx_pic = zeros(size(max_pic2));
    for yy = 1:size(max_pic,1)
         nuer_idx_pic(yy,nuer_idx2(yy,:)) = 1;
    end

    nuer_idx_pic = logical(nuer_idx_pic);
    %subplot(1,2,1)
    %imshowpair(dend_crop,max_pic)
    %title('Max Pic')
    %subplot(1,2,2)
    %imshowpair(dend_crop,nuer_idx_pic(:,2:end))
    %title('Nuer-Idx Tracing')


    % Close nuer-idx for object indentification and easier viewing
    close all
    %figure()
    neur_close = imclose(nuer_idx_pic,strel('rectangle',[11 1]));
    neur_close = neur_close(:,2:end);
    neur_close = bwareaopen(neur_close,3);
    %imshowpair(dend_crop,neur_close)
    %title('Closed Neur-Idx')
    
    % Obtain moving mean (and others) of neuron tracking
    nuer_idx3 = nuer_idx;
    for rr = 1:size(nuer_idx,1)
        for tt = 1:4
            if nuer_idx(rr,tt) == 0
                nuer_idx3(rr,tt) = neur_idx_avg(tt);
            end
        end
    end

    neur_mvmean = movmean(nuer_idx3,meanrange);
    neur_tot_avg = mean(neur_mvmean,1);
    neur_mvvar = movvar(nuer_idx3,20);
    neur_mvmad = movmad(nuer_idx3,20);

    %% Not try to place created 'objects' into correct dendrite
    
    % Get centroids of each object
    neur_close_props = regionprops(neur_close,'all');
    neur_close_cent = regionprops(neur_close,'centroid');
    close_cent = cell2mat(struct2cell(neur_close_cent));
    x_cent = close_cent(1:2:end);
    y_cent = close_cent(2:2:end);
    cent = [x_cent;y_cent]';

    % Create range for binning
    neur_loc_avg(1) = round(mean(neur_idx_avg(1:2)));
    neur_loc_avg(2) = round(mean(neur_idx_avg(2:3)));
    neur_loc_avg(3) = round(mean(neur_idx_avg(3:4)));
    neur_avg_rng(1,:) = [2 neur_loc_avg(1)];
    neur_avg_rng(2,:) = [neur_loc_avg(1) neur_loc_avg(2)];
    neur_avg_rng(3,:) = [neur_loc_avg(2) neur_loc_avg(3)];
    neur_avg_rng(4,:) = [neur_loc_avg(3) size(nuer_idx_pic,2)];

    %% Create a running average for binning
    neur_filt = nuer_idx;
    close all
    img_split = round(size(nuer_idx,1)/split); % Split up the image since x location changes throughout y
    bin_avg = zeros(split,4);
    bin_loc_avg = zeros(split,3);
    bin_avg_rng = zeros(4,2,split);
    y_bin = zeros(split,2);
    binNan = zeros(split,4);

    for run = 1:split
        top = (run-1)*img_split+1:run*img_split;
        topover = find(top>=size(nuer_idx,1));
        top(topover) = [];

        y_bin(run,:) = [top(1)-.99,top(length(top))+.01];

        bin_avg(run,:) = round([sum(nuer_idx(top,1))./nnz(nuer_idx(top,1)), sum(nuer_idx(top,2))./nnz(nuer_idx(top,2)), ...
        sum(nuer_idx(top,3))./nnz(nuer_idx(top,3)), sum(nuer_idx(top,4))./nnz(nuer_idx(top,4))]);

        binNan(run,:) = isnan(bin_avg(run,:));
        for kk = 1:4
           if binNan(run,kk) == 1
               bin_avg(run,kk) = neur_idx_avg(kk);
           end
        end

    % Create range for binning
        bin_loc_avg(run,1) = round(mean(bin_avg(run,1:2)));
        bin_loc_avg(run,2) = round(mean(bin_avg(run,2:3)));
        bin_loc_avg(run,3) = round(mean(bin_avg(run,3:4)));
        bin_avg_rng(1,:,run) = [2 bin_loc_avg(run,1)];
        bin_avg_rng(2,:,run) = [bin_loc_avg(run,1) bin_loc_avg(run,2)];
        bin_avg_rng(3,:,run) = [bin_loc_avg(run,2) bin_loc_avg(run,3)];
        bin_avg_rng(4,:,run) = [bin_loc_avg(run,3) size(nuer_idx_pic,2)];

    end
    
    % Edit moving mean to get out of outliers
    for rr = 1:size(neur_mvmean,1)
        for zone = 1:meansplit
            if rr < y_bin(zone,2) && rr >= y_bin(zone,1)
                xx = floor(y_bin(zone,1)):ceil(y_bin(zone,2));
                xxunder = find(xx<1);
                xx(xxunder) = [];
                xxover = find(xx>size(neur_mvmean,1));
                xx(xxover) = [];
                neur_zone_avg = mean(round(neur_mvmean(xx,:)),1);
                zone_idx = zone;
            end
        end   
        for tt = 1:4
           absdif = abs(neur_mvmean(rr,tt)-neur_zone_avg(tt));
           if absdif > 7
               neur_mvmean(rr,tt) = neur_zone_avg(tt);
           end
        end
    end

    % Place object in bins based on x value of centroid
    step = 0;
    binIdx = zeros(length(x_cent),4);
    for ii = 1:length(x_cent)
        for zone = 1:split
            if y_cent(ii) < y_bin(zone,2) && y_cent(ii) >= y_bin(zone,1)
                zone_idx = zone;
            end
        end

        for jj = 1:4
            if x_cent(ii) < bin_avg_rng(jj,2,zone_idx) && x_cent(ii) >= bin_avg_rng(jj,1,zone_idx)
                binIdx(ii,jj) = ii;
            end
        end 

    end

    dend1idx = binIdx(:,1);
    dend1idx(dend1idx==0) = [];
    dend2idx = binIdx(:,2);
    dend2idx(dend2idx==0) = [];
    dend3idx = binIdx(:,3);
    dend3idx(dend3idx==0) = [];
    dend4idx = binIdx(:,4);
    dend4idx(dend4idx==0) = [];

    % Create images with binned objects
    dend1close = neur_close;
    dend2close = neur_close;
    dend3close = neur_close;
    dend4close = neur_close;

    for uu = 1:length(x_cent)
        if uu ~= dend1idx
           dend1close(neur_close_props(uu).PixelIdxList) = 0;
        end
        if uu ~= dend2idx
            dend2close(neur_close_props(uu).PixelIdxList) = 0;
        end
        if uu ~= dend3idx
            dend3close(neur_close_props(uu).PixelIdxList) = 0;
        end
        if uu ~= dend4idx
            dend4close(neur_close_props(uu).PixelIdxList) = 0;
        end
    end
    %figure(1)
    %subplot(1,4,1)
    %imshowpair(dend_crop,dend1close)
    %title('Dendrite 1')
    %subplot(1,4,2)
    %imshowpair(dend_crop,dend2close)
    %title('Dendrite 2')
    %subplot(1,4,3)
    %imshowpair(dend_crop,dend3close)
    %title('Dendrite 3')
    %subplot(1,4,4)
    %imshowpair(dend_crop,dend4close)
    %title('Dendrite 4')
    %sgtitle('Original Binning - Closed')

    %% Attempt skeletonizing for ranges
    close all
    dend1skel = bwskel(dend1close,'MinBranchLength',2);
    dend2skel = bwskel(dend2close,'MinBranchLength',2);
    dend3skel = bwskel(dend3close,'MinBranchLength',2);
    dend4skel = bwskel(dend4close,'MinBranchLength',2);
    dend1_lbl = bwlabel(dend1skel);
    dend2_lbl = bwlabel(dend2skel);
    dend3_lbl = bwlabel(dend3skel);
    dend4_lbl = bwlabel(dend4skel);
    dend1cent = regionprops(dend1skel,'centroid');
    dend1cent = cell2mat(struct2cell(dend1cent));
    dend1cent = dend1cent(1:2:end);
    dend2cent = regionprops(dend2skel,'centroid');
    dend2cent = cell2mat(struct2cell(dend2cent));
    dend2cent = dend2cent(1:2:end);
    dend3cent = regionprops(dend3skel,'centroid');
    dend3cent = cell2mat(struct2cell(dend3cent));
    dend3cent = dend3cent(1:2:end);
    dend4cent = regionprops(dend4skel,'centroid');
    dend4cent = cell2mat(struct2cell(dend4cent));
    dend4cent = dend4cent(1:2:end);

    %{ 
    subplot(1,4,1)
    %imshowpair(dend_crop,dend1skel)
    title('Dendrite 1')
    subplot(1,4,2)
    %imshowpair(dend_crop,dend2skel)
    title('Dendrite 2')
    subplot(1,4,3)
    %imshowpair(dend_crop,dend3skel)
    title('Dendrite 3')
    subplot(1,4,4)
    %imshowpair(dend_crop,dend4skel)
    title('Dendrite 4')
    sgtitle('Original Binning - Skletonized')
    %} 
    
    
    %% Remove pixels farther from the average
    % close all
    % Dendrite 1
    obj2mov1up = zeros(size(dend_crop));
    for yy = 2:size(dend_crop,1)-1
        dend1_temp = dend1skel(yy,:);
        tempidx = find(dend1_temp);
        if length(tempidx) > 1
            for zone = 1:split
                if yy < y_bin(zone,2) && yy >= y_bin(zone,1)
                    zone_idx = zone;
                end
            end
            tempobj = dend1_lbl(yy,tempidx);
            tempdist = zeros(length(tempidx),1);

            for zz = 1:length(tempidx)
                tempdist(zz) = abs(tempidx(zz)-bin_avg(zone_idx,1));
            end
                [mindist, mindist_idx] = min(tempdist);
                real_idx = tempidx(mindist_idx);
                real_obj = tempobj(mindist_idx);

            for i = 1:length(tempidx)
                tempup = zeros(1,3);
                tempdown = zeros(1,3);
                tempup(1,:) = dend1_lbl(yy-1,tempidx(i)-1:tempidx(i)+1);
                tempdown(1,:) = dend1_lbl(yy+1,tempidx(i)-1:tempidx(i)+1);
                tempup(tempup<=0) = [];
                tempdown(tempdown<=0) = [];
                if isempty(tempup) ~= 1 || isempty(tempdown) ~= 1
                    if isempty(tempup) == 1
                        tempup = 0;
                    end
                    if isempty(tempdown) == 1
                        tempdown = 0;
                    end
                    if tempobj(i) == tempup(1) || tempobj(i) == tempdown(1)
                        temp_cent = dend1cent(tempobj(i));
                        if temp_cent < neur_avg_rng(1,1) || temp_cent > neur_avg_rng(1,2)
                            if tempidx(i) ~= real_idx && tempidx(i) > real_idx && tempobj(i) ~= real_obj
                                obj2mov1up(yy,tempidx(i)) = dend1_lbl(yy,tempidx(i));
                                dend1skel(yy,tempidx(i)) = 0;
                            end
                        end

                    end
                elseif tempidx(i) ~= real_idx && tempidx(i) > real_idx && tempobj(i) ~= real_obj
                        obj2mov1up(yy,tempidx(i)) = dend1_lbl(yy,tempidx(i));
                        dend1skel(yy,tempidx(i)) = 0;

                end            

            end
        end
    end
    [obj1uprow, obj1upcol] = find(obj2mov1up);
    for i = 1:length(obj1uprow)
        dend2skel(obj1uprow(i), obj1upcol(i)) = 1;
    end

    % Dendrite 2
    obj2mov2up = zeros(size(dend_crop));
    obj2mov2dn = zeros(size(dend_crop));
    for yy = 2:size(dend_crop,1)-1
        dend2_temp = dend2skel(yy,:);
        tempidx = find(dend2_temp);
        if length(tempidx) > 1
            for zone = 1:split
                if yy < y_bin(zone,2) && yy >= y_bin(zone,1)
                    zone_idx = zone;
                end
            end
            tempobj = dend2_lbl(yy,tempidx);
            tempdist = zeros(length(tempidx),1);

            for zz = 1:length(tempidx)
                tempdist(zz) = abs(tempidx(zz)-bin_avg(zone_idx,2));
            end
                [mindist, mindist_idx] = min(tempdist);
                real_idx = tempidx(mindist_idx);
                real_obj = tempobj(mindist_idx);

            for i = 1:length(tempidx)
                tempup = zeros(1,3);
                tempdown = zeros(1,3);
                tempup(1,:) = dend2_lbl(yy-1,tempidx(i)-1:tempidx(i)+1);
                tempdown(1,:) = dend2_lbl(yy+1,tempidx(i)-1:tempidx(i)+1);
                tempup(tempup<=0) = [];
                tempdown(tempdown<=0) = [];
                if isempty(tempup) ~= 1 || isempty(tempdown) ~= 1
                    if isempty(tempup) == 1
                        tempup = 0;
                    end
                    if isempty(tempdown) == 1
                        tempdown = 0;
                    end
                    if tempobj(i) == tempup(1) || tempobj(i) == tempdown(1)
                        temp_cent = dend2cent(tempobj(i));
                        if temp_cent < neur_avg_rng(2,1) || temp_cent > neur_avg_rng(2,2)
                            if tempidx(i) ~= real_idx && tempidx(i) > real_idx && tempobj(i) ~= real_obj
                                obj2mov2up(yy,tempidx(i)) = dend2_lbl(yy,tempidx(i));
                                dend2skel(yy,tempidx(i)) = 0;
                            elseif tempidx(i) ~= real_idx && tempidx(i) < real_idx && tempobj(i) ~= real_obj
                                obj2mov2dn(yy,tempidx(i)) = dend2_lbl(yy,tempidx(i));
                                dend2skel(yy,tempidx(i)) = 0;
                            end
                        end

                    end
                elseif tempidx(i) ~= real_idx && tempidx(i) > real_idx && tempobj(i) ~= real_obj
                    obj2mov2up(yy,tempidx(i)) = dend2_lbl(yy,tempidx(i));
                    dend2skel(yy,tempidx(i)) = 0;
                elseif tempidx(i) ~= real_idx && tempidx(i) < real_idx && tempobj(i) ~= real_obj
                    obj2mov2dn(yy,tempidx(i)) = dend2_lbl(yy,tempidx(i));
                    dend2skel(yy,tempidx(i)) = 0;
                end            

            end
        end
    end
    [obj2uprow, obj2upcol] = find(obj2mov2up);
    [obj2dnrow, obj2dncol] = find(obj2mov2dn);
    for i = 1:length(obj2uprow)
        dend3skel(obj2uprow(i), obj2upcol(i)) = 1;
    end
    for i = 1:length(obj2dnrow)
        dend1skel(obj2dnrow(i), obj2dncol(i)) = 1;
    end

    % Dendrite 3
    obj2mov3up = zeros(size(dend_crop));
    obj2mov3dn = zeros(size(dend_crop));
    for yy = 2:size(dend_crop,1)-1
        dend3_temp = dend3skel(yy,:);
        tempidx = find(dend3_temp);
        if length(tempidx) > 1
            for zone = 1:split
                if yy < y_bin(zone,2) && yy >= y_bin(zone,1)
                    zone_idx = zone;
                end
            end
            tempobj = dend3_lbl(yy,tempidx);
            tempdist = zeros(length(tempidx),1);
            for zz = 1:length(tempidx)
                tempdist(zz) = abs(tempidx(zz)-bin_avg(zone_idx,3));
            end
                [mindist, mindist_idx] = min(tempdist);
                real_idx = tempidx(mindist_idx);
                real_obj = tempobj(mindist_idx);
            for i = 1:length(tempidx)
                tempup = zeros(1,3);
                tempdown = zeros(1,3);
                tempup(1,:) = dend3_lbl(yy-1,tempidx(i)-1:tempidx(i)+1);
                tempdown(1,:) = dend3_lbl(yy+1,tempidx(i)-1:tempidx(i)+1);
                tempup(tempup<=0) = [];
                tempdown(tempdown<=0) = [];
                if isempty(tempup) ~= 1 || isempty(tempdown) ~= 1
                    if isempty(tempup) == 1
                        tempup = 0;
                    end
                    if isempty(tempdown) == 1
                        tempdown = 0;
                    end
                    if tempobj(i) == tempup(1) || tempobj(i) == tempdown(1)
                        temp_cent = dend3cent(tempobj(i));
                        if temp_cent < neur_avg_rng(3,1) || temp_cent > neur_avg_rng(3,2)
                            if tempidx(i) ~= real_idx && tempidx(i) > real_idx && tempobj(i) ~= real_obj
                                obj2mov3up(yy,tempidx(i)) = dend3_lbl(yy,tempidx(i));
                                dend3skel(yy,tempidx(i)) = 0;
                            elseif tempidx(i) ~= real_idx && tempidx(i) < real_idx && tempobj(i) ~= real_obj
                                obj2mov3dn(yy,tempidx(i)) = dend3_lbl(yy,tempidx(i));
                                dend3skel(yy,tempidx(i)) = 0;
                            end
                        end

                    end
                elseif tempidx(i) ~= real_idx && tempidx(i) > real_idx && tempobj(i) ~= real_obj
                    obj2mov3up(yy,tempidx(i)) = dend3_lbl(yy,tempidx(i));
                    dend3skel(yy,tempidx(i)) = 0;
                elseif tempidx(i) ~= real_idx && tempidx(i) < real_idx && tempobj(i) ~= real_obj
                    obj2mov3dn(yy,tempidx(i)) = dend3_lbl(yy,tempidx(i));
                    dend3skel(yy,tempidx(i)) = 0;
                end            

            end
        end
    end
    [obj3uprow, obj3upcol] = find(obj2mov3up);
    [obj3dnrow, obj3dncol] = find(obj2mov3dn);
    for i = 1:length(obj3uprow)
        dend4skel(obj3uprow(i), obj3upcol(i)) = 1;
    end
    for i = 1:length(obj3dnrow)
        dend2skel(obj3dnrow(i), obj3dncol(i)) = 1;
    end

    % Dendrite 4
    obj2mov4up = zeros(size(dend_crop));
    obj2mov4dn = zeros(size(dend_crop));
    for yy = 2:size(dend_crop,1)-1
        dend4_temp = dend4skel(yy,:);
        tempidx = find(dend4_temp);
        if length(tempidx) > 1
            for zone = 1:split
                if yy < y_bin(zone,2) && yy >= y_bin(zone,1)
                    zone_idx = zone;
                end
            end
            tempobj = dend4_lbl(yy,tempidx);
            tempdist = zeros(length(tempidx),1);
            for zz = 1:length(tempidx)
                tempdist(zz) = abs(tempidx(zz)-bin_avg(zone_idx,4));
            end
                [mindist, mindist_idx] = min(tempdist);
                real_idx = tempidx(mindist_idx);
                real_obj = tempobj(mindist_idx);
            for i = 1:length(tempidx)  
                tempup = zeros(1,3);
                tempdown = zeros(1,3);
                tempup(1,:) = dend4_lbl(yy-1,tempidx(i)-1:tempidx(i)+1);
                tempdown(1,:) = dend4_lbl(yy+1,tempidx(i)-1:tempidx(i)+1);
                tempup(tempup<=0) = [];
                tempdown(tempdown<=0) = [];
                if isempty(tempup) ~= 1 || isempty(tempdown) ~= 1
                    if isempty(tempup) == 1
                        tempup = 0;
                    end
                    if isempty(tempdown) == 1
                        tempdown = 0;
                    end
                    if tempobj(i) == tempup(1) || tempobj(i) == tempdown(1)
                        temp_cent = dend4cent(tempobj(i));
                        if temp_cent < neur_avg_rng(4,1) || temp_cent > neur_avg_rng(4,2)
                            if tempidx(i) ~= real_idx && tempidx(i) < real_idx && tempobj(i) ~= real_obj
                                obj2mov4dn(yy,tempidx(i)) = dend4_lbl(yy,tempidx(i));
                                dend4skel(yy,tempidx(i)) = 0;
                            end
                        end
                    end
                elseif tempidx(i) ~= real_idx && tempidx(i) < real_idx && tempobj(i) ~= real_obj
                    obj2mov4dn(yy,tempidx(i)) = dend4_lbl(yy,tempidx(i));
                    dend4skel(yy,tempidx(i)) = 0;
                end            

            end
        end
    end
    [obj4dnrow, obj4dncol] = find(obj2mov4dn);
    for i = 1:length(obj4dnrow)
        dend3skel(obj4dnrow(i), obj4dncol(i)) = 1;
    end

    %{
    % View skeletonized version
    figure(2)
    subplot(1,4,1)
    %imshowpair(dend_crop,dend1skel)
    title('Dendrite 1')
    subplot(1,4,2)
    %imshowpair(dend_crop,dend2skel)
    title('Dendrite 2')
    subplot(1,4,3)
    %imshowpair(dend_crop,dend3skel)
    title('Dendrite 3')
    subplot(1,4,4)
    %imshowpair(dend_crop,dend4skel)
    title('Dendrite 4')
    sgtitle('After Moving Points')
    %} 
    
    %% Run through original binning to get rid of duplicates
    close all
    dend1_lbl = bwlabel(dend1skel);
    dend2_lbl = bwlabel(dend2skel);
    dend3_lbl = bwlabel(dend3skel);
    dend4_lbl = bwlabel(dend4skel);

    dend1area = regionprops(dend1skel,'area');
    dend1area = [dend1area.Area];
    [~,dend1max] = max(dend1area);
    dend1cent = regionprops(dend1skel,'centroid');
    dend1cent = cell2mat(struct2cell(dend1cent));
    dend1cent = dend1cent(1:2:end);

    dend2area = regionprops(dend2skel,'area');
    dend2area = [dend2area.Area];
    [~,dend2max] = max(dend2area);
    dend2cent = regionprops(dend2skel,'centroid');
    dend2cent = cell2mat(struct2cell(dend2cent));
    dend2cent = dend2cent(1:2:end);

    dend3area = regionprops(dend3skel,'area');
    dend3area = [dend3area.Area];
    [~,dend3max] = max(dend3area);
    dend3cent = regionprops(dend3skel,'centroid');
    dend3cent = cell2mat(struct2cell(dend3cent));
    dend3cent = dend3cent(1:2:end);

    dend4area = regionprops(dend4skel,'area');
    dend4area = [dend4area.Area];
    [~,dend4max] = max(dend4area);
    dend4cent = regionprops(dend4skel,'centroid');
    dend4cent = cell2mat(struct2cell(dend4cent));
    dend4cent = dend4cent(1:2:end);

    % SinglePointBin Function Run
    dend1skel = SinglePointBin(dend1skel,dend1_lbl,dend1cent,dend1max,1,neur_mvmean,neur_idx_avg);
    dend2skel = SinglePointBin(dend2skel,dend2_lbl,dend2cent,dend2max,2,neur_mvmean,neur_idx_avg);
    dend3skel = SinglePointBin(dend3skel,dend3_lbl,dend3cent,dend3max,3,neur_mvmean,neur_idx_avg);
    dend4skel = SinglePointBin(dend4skel,dend4_lbl,dend4cent,dend4max,4,neur_mvmean,neur_idx_avg);

    %{ 
    % View skeletonized version
    subplot(1,4,1)
    %imshowpair(dend_crop,dend1skel)
    title('Dendrite 1')
    subplot(1,4,2)
    %imshowpair(dend_crop,dend2skel)
    title('Dendrite 2')
    subplot(1,4,3)
    %imshowpair(dend_crop,dend3skel)
    title('Dendrite 3')
    subplot(1,4,4)
    %imshowpair(dend_crop,dend4skel)
    title('Dendrite 4')
    sgtitle('Binning After SinglePointBin')
    %} 
    
    %% Using moving MAD to remove weird jumps
    dend1skel2 = dend1skel;
    dend2skel2 = dend2skel;
    dend3skel2 = dend3skel;
    dend4skel2 = dend4skel;
    dend1skel_idx = zeros(size(dend1skel2,1),1);
    dend2skel_idx = zeros(size(dend2skel2,1),1);
    dend3skel_idx = zeros(size(dend3skel2,1),1);
    dend4skel_idx = zeros(size(dend4skel2,1),1);
    dend1id = zeros(size(dend1skel2,1),1);
    dend2id = zeros(size(dend2skel2,1),1);
    dend3id = zeros(size(dend3skel2,1),1);
    dend4id = zeros(size(dend4skel2,1),1);
    for rr = 1:size(dend4skel2,1)
        dend1temp = find(dend1skel2(rr,:));
        if isempty(dend1temp) == 1
            dend1skel_idx(rr) = neur_mvmean(rr,1);
            dend1id(rr) = 0;
        elseif length(dend1temp) ==1
            dend1skel_idx(rr) = find(dend1skel2(rr,:));
            dend1id(rr) = dend1temp;
        elseif length(dend1temp) > 1
            dend1skel_idx(rr) = mean(find(dend1skel2(rr,:)));
            dend1id(rr) = mean(dend1temp);
        end
        dend2temp = find(dend2skel2(rr,:));
        if isempty(dend2temp) == 1
            dend2skel_idx(rr) = neur_mvmean(rr,2);
            dend2id(rr) = 0;
        elseif length(dend2temp) ==1
            dend2skel_idx(rr) = find(dend2skel2(rr,:));
            dend2id(rr) = dend2temp;
        elseif length(dend2temp) > 1
            dend2skel_idx(rr) = mean(find(dend2skel2(rr,:)));
            dend2id(rr) = mean(dend2temp);
        end
        dend3temp = find(dend3skel2(rr,:));
        if isempty(dend3temp) == 1
            dend3skel_idx(rr) = neur_mvmean(rr,3);
            dend3id(rr) = 0;
        elseif length(dend3temp) ==1
            dend3skel_idx(rr) = find(dend3skel2(rr,:));
            dend3id(rr) = dend3temp;
        elseif length(dend3temp) > 1
            dend3skel_idx(rr) = mean(find(dend3skel2(rr,:)));
            dend3id(rr) = mean(dend3temp);
        end
        dend4temp = find(dend4skel2(rr,:));
        if isempty(dend4temp) == 1
            dend4skel_idx(rr) = neur_mvmean(rr,4);
            dend4id(rr) = 0;
        elseif length(dend4temp) ==1
            dend4skel_idx(rr) = find(dend4skel2(rr,:));
            dend4id(rr) = dend4temp;
        elseif length(dend4temp) > 1
            dend4skel_idx(rr) = mean(find(dend4skel2(rr,:)));
            dend4id(rr) = mean(dend4temp);
        end

            dend1skel2(rr,round(dend1skel_idx(rr))) = 1;
            dend2skel2(rr,round(dend2skel_idx(rr))) = 1;
            dend3skel2(rr,round(dend3skel_idx(rr))) = 1;
            dend4skel2(rr,round(dend4skel_idx(rr))) = 1;
    end
    % Make images and plot!
    dend1 = zeros(size(dend1skel));
    dend2 = zeros(size(dend2skel));
    dend3 = zeros(size(dend3skel));
    dend4 = zeros(size(dend4skel));
    dend1id(dend1id<=1) = 1;
    dend2id(dend2id<=1) = 1;
    dend3id(dend3id<=1) = 1;
    dend4id(dend4id<=1) = 1;

    for rr = 1:length(dend1id)
        dend1(rr,dend1id(rr)) = 1;
        dend2(rr,dend2id(rr)) = 1;
        dend3(rr,dend3id(rr)) = 1;
        dend4(rr,dend4id(rr)) = 1;    
    end
    dend1 = logical(dend1(:,2:end));
    dend2 = logical(dend2(:,2:end));
    dend3 = logical(dend3(:,2:end));
    dend4 = logical(dend4(:,2:end));
    
    %{ 
    % View skeletonized version
    subplot(1,4,1)
    %imshowpair(dend_crop,dend1skel2)
    title('Dendrite 1')
    subplot(1,4,2)
    %imshowpair(dend_crop,dend2skel2)
    title('Dendrite 2')
    subplot(1,4,3)
    %imshowpair(dend_crop,dend3skel2)
    title('Dendrite 3')
    subplot(1,4,4)
    %imshowpair(dend_crop,dend4skel2)
    title('Dendrite 4')
    sgtitle('Binning After Removing Weird Jumps','fontsize',26)
    %} 

    %% Try to smooth data (using smooth data) and dend skel
    close all
    [dend1smooth, dend1smoothidx] = NeurSmoothData(dend1id,dend1skel2);
    [dend2smooth, dend2smoothidx] = NeurSmoothData(dend2id,dend2skel2);
    [dend3smooth, dend3smoothidx, dend3b4nan, dend3nan] = NeurSmoothData(dend3id,dend3skel2);
    [dend4smooth, dend4smoothidx, dend4b4nan, dend4nan] = NeurSmoothData(dend4id,dend4skel2);
    
    %{ 
    subplot(2,4,1)
    %imshowpair(dend_crop,dend1smooth)
    title('Dendrite 1')
    subplot(2,4,2)
    %imshowpair(dend_crop,dend2smooth)
    title('Dendrite 2')
    subplot(2,4,3)
    %imshowpair(dend_crop,dend3smooth)
    title('Dendrite 3')
    subplot(2,4,4)
    %imshowpair(dend_crop,dend4smooth)
    title('Dendrite 4')

    subplot(2,4,5)
    %imshowpair(dend_crop,dend1)
    title('Dendrite 1')
    subplot(2,4,6)
    %imshowpair(dend_crop,dend2)
    title('Dendrite 2')
    subplot(2,4,7)
    %imshowpair(dend_crop,dend3)
    title('Dendrite 3')
    subplot(2,4,8)
    %imshowpair(dend_crop,dend4)
    title('Dendrite 4')
    sgtitle('After Data Smoothing!');
    %} 

    %% Now interpolate for each section
    dend1_idx = zeros(size(dend_crop,1),1);
    dend2_idx = zeros(size(dend_crop,1),1);
    dend3_idx = zeros(size(dend_crop,1),1);
    dend4_idx = zeros(size(dend_crop,1),1);

    % Index each dendrite
    for xx = 1:size(dend_crop,1)
        if isempty(find(dend1smooth(xx,:),1)) == 1
            dend1_idx(xx) = 0;
        else
            dend1_idx(xx) = find(dend1smooth(xx,:));
        end
        if isempty(find(dend2smooth(xx,:),1)) == 1
            dend2_idx(xx) = 0;
        else
            dend2_idx(xx) = find(dend2smooth(xx,:));
        end    
        if isempty(find(dend3smooth(xx,:),1)) == 1
            dend3_idx(xx) = 0;
        else
            dend3_idx(xx) = find(dend3smooth(xx,:));
        end   
        if isempty(find(dend4smooth(xx,:),1)) == 1
            dend4_idx(xx) = 0;
        else
            dend4_idx(xx) = find(dend4smooth(xx,:));
        end    

    end
    idx1 = find(dend1_idx);
    idx2 = find(dend2_idx);
    idx3 = find(dend3_idx);
    idx4 = find(dend4_idx);

    % Interpolate
    neur1 = dend1_idx;
    nz_neur1 = find(neur1~=0);
    nz_input1 = neur1(nz_neur1);
    neur1 = interp1(nz_neur1,nz_input1,1:length(neur1),'linear','extrap');
    neur1(isnan(neur1)) = 0;
    nuer(:,1) = neur1;

    neur2 = dend2_idx;
    nz_neur2 = find(neur2~=0);
    nz_input2 = neur2(nz_neur2);
    neur2 = interp1(nz_neur2,nz_input2,1:length(neur2),'linear','extrap');
    neur2(isnan(neur2)) = 0;
    nuer(:,2) = neur2;

    neur3 = dend3_idx;
    nz_neur3 = find(neur3~=0);
    nz_input3 = neur3(nz_neur3);
    neur3 = interp1(nz_neur3,nz_input3,1:length(neur3),'linear','extrap');
    neur3(isnan(neur3)) = 0;
    nuer(:,3) = neur3;

    neur4 = dend4_idx;
    nz_neur4 = find(neur4~=0);
    nz_input4 = neur4(nz_neur4);
    if idx4(1) > length(neur4)/2 
        neur4 = interp1(nz_neur4,nz_input4,1:length(neur4),'linear',dend4smoothidx(idx4(1)));
    elseif idx4(end) < length(neur4)/2
        neur4 = interp1(nz_neur4,nz_input4,1:length(neur4),'linear',dend4smoothidx(idx4(end)));
    else
        neur4 = interp1(nz_neur4,nz_input4,1:length(neur4),'linear','extrap');
    end
    neur4(isnan(neur4)) = 0;
    nuer(:,4) = neur4;

    nuer = round(nuer);
    
    %% Create individual dendrite range 
    nuer_range = zeros(size(nuer_idx,1),2,4);
    obsRng = 10;
    for zz = 1:size(nuer_idx,1)
        for mm = 1:4
            nuer_range(zz,:,mm) = [nuer(zz,mm)-obsRng nuer(zz,mm)+obsRng];
        end
    end
    nuer_range(nuer_range<=0) = 1;
    nuer_range(nuer_range>size(dend_crop,2)) = size(dend_crop,2);

    %% threshold within each neuron range and creat full 'mask'
    close all
    nuer1mask = zeros(size(dend_crop));
    nuer2mask = zeros(size(dend_crop));
    nuer3mask = zeros(size(dend_crop));
    nuer4mask = zeros(size(dend_crop));

    for ll = 1:size(dend_crop,1)
        nuer1mask(ll,nuer_range(ll,1,1):nuer_range(ll,2,1)) = 1;
        nuer2mask(ll,nuer_range(ll,1,2):nuer_range(ll,2,2)) = 1;
        nuer3mask(ll,nuer_range(ll,1,3):nuer_range(ll,2,3)) = 1;
        nuer4mask(ll,nuer_range(ll,1,4):nuer_range(ll,2,4)) = 1;
    end

    nuermask = nuer1mask | nuer2mask | nuer3mask | nuer4mask;

    %imshowpair(dend_crop,nuermask)
    
    
    %% Match to reference image and then crop dend_crop to make sure range is used
    dendcrop_match = imhistmatch(dend_crop,dend_crop);
    %imshowpair(dend_crop,dendcrop_match,'montage')

    %% Find first point to avoid erros from rotating
    firstPt = zeros(size(dend_crop,2),1);
    for ii = 1:size(dend_crop,2)
        firstPt(ii) = find(dend_crop(:,ii),1,'first');
    end
    firstPt = max(firstPt);


    dend_crop2 = imcrop(dendcrop_match,[1,firstPt,size(dend_crop,2),size(dend_crop,1)-firstPt]);
    nuer1mask2 = imcrop(nuer1mask,[1,firstPt,size(dend_crop,2),size(dend_crop,1)-firstPt]);
    nuer2mask2 = imcrop(nuer2mask,[1,firstPt,size(dend_crop,2),size(dend_crop,1)-firstPt]);
    nuer3mask2 = imcrop(nuer3mask,[1,firstPt,size(dend_crop,2),size(dend_crop,1)-firstPt]);
    nuer4mask2 = imcrop(nuer4mask,[1,firstPt,size(dend_crop,2),size(dend_crop,1)-firstPt]);
    nuermask2 = imcrop(nuermask,[1,firstPt,size(dend_crop,2),size(dend_crop,1)-firstPt]);


    %% Mask out dendrite regions
    close all
    dend1img = bsxfun(@times, dend_crop2, cast(nuer1mask2, 'like', dend_crop2)); % mask out image
    dend2img = bsxfun(@times, dend_crop2, cast(nuer2mask2, 'like', dend_crop2));
    dend3img = bsxfun(@times, dend_crop2, cast(nuer3mask2, 'like', dend_crop2));
    dend4img = bsxfun(@times, dend_crop2, cast(nuer4mask2, 'like', dend_crop2));
    dendimg = bsxfun(@times, dend_crop2, cast(nuermask2, 'like', dend_crop2));

    %% NEW CODE (11/10/2020) TO RESHAPE DEND IMAGES
    % d1 = dend1img;
    % d1 = d1';
    % d1(d1==0) = [];
    % d1 = imresize(reshape(d1,11,size(dend1img,1))',4);
    d1f = imflatfield(dend1img, [50 size(dend4img,2)]);
    % d2 = dend2img;
    % d2 = d2';
    % d2(d2==0) = [];
    % d2 = imresize(reshape(d2,11,size(dend2img,1))',4);
    d2f = imflatfield(dend2img,[50 size(dend4img,2)]);
    % d3 = dend3img;
    % d3 = d3';
    % d3(d3==0) = [];
    % d3 = imresize(reshape(d3,11,size(dend3img,1))',4);
    d3f = imflatfield(dend3img,[50 size(dend4img,2)]);
    % d4 = dend4img;
    % d4 = d4';
    % d4(d4==0) = [];
    % d4 = imresize(reshape(d4,11,size(dend4img,1))',4);
    d4f = imflatfield(dend4img,[50 size(dend4img,2)]);
    %% Lets get them blebs and breaks
    % Attempt to analyze each section using islocalmax
    close all
    dend1width = zeros(size(dend1img));
    dend2width = zeros(size(dend2img));
    dend3width = zeros(size(dend3img));
    dend4width = zeros(size(dend4img));
    dendthresh = double(prctile(nonzeros(dendimg),25));
    dendthreshbleb = double(prctile(nonzeros(dendimg),95));
    dendthreshtot = double(prctile(nonzeros(dend_crop2),25));
    blebPts = zeros(size(dend1img,1),4);
    [blebPts1, width1, oneidx1, dendidx1, S1, dendbw1, dendtest, chg1, mintempS1, dendthresh1] = blebIntensityCount(dendimg,dend1img,1,dendthreshtot);
    [blebPts2, width2, oneidx2, dendidx2, S2, dendbw2, dendtest2, chg2, mintempS2, dendthresh2] = blebIntensityCount(dendimg,dend2img,2,dendthreshtot);
    [blebPts3, width3, oneidx3, dendidx3, S3, dendbw3, dendtest3, chg3, mintempS3, dendthresh3] = blebIntensityCount(dendimg,dend3img,3,dendthreshtot);
    [blebPts4, width4, oneidx4, dendidx4, S4, dendbw4, dendtest4, chg4, mintempS4, dendthresh4]  = blebIntensityCount(dendimg,dend4img,4,dendthreshtot);
    set(gcf,'Position',[1000,10,1000,1000])

    %% Get row intensity average for each dendrite 
    [rowAvg1, rowAvgMean1, rowAvgSTD1] = getRowIntAvg(dend1img,dendbw1);
    [rowAvg2, rowAvgMean2, rowAvgSTD2] = getRowIntAvg(dend2img,dendbw2);
    [rowAvg3, rowAvgMean3, rowAvgSTD3] = getRowIntAvg(dend3img,dendbw3);
    [rowAvg4, rowAvgMean4, rowAvgSTD4] = getRowIntAvg(dend4img,dendbw4);
    rowImgMed = median([rowAvgMean1,rowAvgMean2,rowAvgMean3,rowAvgMean4]);


    %% Average dendbw and intensity to get blebs
    close all
    [blebLoc1, blebCount1, maxIntProm1, blebInt1, brkcount1, brkprct1, xbreakLoc1, findbreakLoc1, thresh1, rowIntAvg1, blebOut1] = findBlebsAndBreaks(dend1img,dendbw1,dendidx1,1,img_num);
    [blebLoc2, blebCount2, maxIntProm2, blebInt2, brkcount2, brkprct2, xbreakLoc2, findbreakLoc2, thresh2, rowIntAvg2, blebOut2] = findBlebsAndBreaks(dend2img,dendbw2,dendidx2,2,img_num);
    [blebLoc3, blebCount3, maxIntProm3, blebInt3, brkcount3, brkprct3, xbreakLoc3, findbreakLoc3, thresh3, rowIntAvg3, blebOut3] = findBlebsAndBreaks(dend3img,dendbw3,dendidx3,3,img_num);
    [blebLoc4, blebCount4, maxIntProm4, blebInt4, brkcount4, brkprct4, xbreakLoc4, findbreakLoc4, thresh4, rowIntAvg4, blebOut4] = findBlebsAndBreaks(dend4img,dendbw4,dendidx4,4,img_num);
    set(gcf,'Position',[1000,10,1000,1000])
    brk = mean([brkprct1, brkprct2,brkprct3, brkprct4]); % Combine breaks
    breakpct = [brkprct1;brkprct2;brkprct3;brkprct4];

    %% Show dendrite ranges
    %{
    close all
    subplot(1,4,1)
    %imshowpair(dend_crop,nuer1mask);
    title('Dendrite 1')

    subplot(1,4,2)
    %imshowpair(dend_crop,nuer2mask);
    title('Dendrite 2')

    subplot(1,4,3)
    %imshowpair(dend_crop,nuer3mask);
    title('Dendrite 3')

    subplot(1,4,4)
    %imshowpair(dend_crop,nuer4mask);
    title('Dendrite 4')
    %} 
    
   %% Identify which blebs on which dendrites 
   
   
   % Determine if objects are in which dendrite region
    overLap1 = nuer1mask & blebMask;
    overLap2 = nuer2mask & blebMask;
    overLap3 = nuer3mask & blebMask;
    overLap4 = nuer4mask & blebMask; 
    
    [blebRow1,blebCol1] = find(overLap1); % List of the coordinates of all found blebs 
    centers1 = [blebCol1 blebRow1]; 
    [blebRow2,blebCol2] = find(overLap2); 
    centers2 = [blebRow2 blebCol2];
    [blebRow3,blebCol3] = find(overLap3);
    centers3 = [blebRow3 blebCol3];
    [blebRow4,blebCol4] = find(overLap4);
    centers4 = [blebRow4 blebCol4]; 
    
    finalBlebCount = length(blebRow1) + length(blebRow2) + length(blebRow3) + length(blebRow4); 
    
    
    %% Show dendrite ranges
    % Put dendrite health in title of each sub image 
    
    
    imNumber = string(img_num); %% Changes img number into a string so I can save file name 
    breakpercent1 = string(brkprct1); 
    breakpercent2 = string(brkprct2); 
    breakpercent3 = string(brkprct3); 
    breakpercent4 = string(brkprct4); 
    
    close all
    subplot(1,4,1)
    imshow(dend_crop);
    hold on 
    plot(blebCol1,blebRow1,'o', 'LineWidth', .5, 'MarkerSize', 11, 'MarkerEdgeColor','g'); 
    plot(blebCol1,blebRow1,'.', 'LineWidth', .5, 'MarkerSize', 1, 'MarkerEdgeColor','g'); 
    title('Dendrite 1 (health: ' + breakpercent1 + '%)') 
    
    
    subplot(1,4,2)
    imshow(dend_crop);
    hold on 
    plot(blebCol2,blebRow2,'o', 'LineWidth', .5, 'MarkerSize', 11, 'MarkerEdgeColor','g');
    plot(blebCol2,blebRow2,'.', 'LineWidth', .5, 'MarkerSize', 1, 'MarkerEdgeColor','g'); 
    title('Dendrite 2 (health: ' + breakpercent2 + '%)') 

    subplot(1,4,3)
    imshow(dend_crop);
    hold on 
    plot(blebCol3,blebRow3,'o', 'LineWidth', .5, 'MarkerSize', 11, 'MarkerEdgeColor','g');
    plot(blebCol3,blebRow3,'.', 'LineWidth', .5, 'MarkerSize', 1, 'MarkerEdgeColor','g'); 
    title('Dendrite 3 (health: ' + breakpercent3 + '%)') 

    subplot(1,4,4)
    imshow(dend_crop);
    hold on 
    plot(blebCol4,blebRow4,'o', 'LineWidth', .5, 'MarkerSize', 11, 'MarkerEdgeColor','g');
    plot(blebCol4,blebRow4,'.', 'LineWidth', .5, 'MarkerSize', 1, 'MarkerEdgeColor','g'); 
    title('Dendrite 4 (health: ' + breakpercent4 + '%)') 
    
    set(gcf, 'Position', [10 10 900 600]); 
    saveas(gcf, strcat('C:\Users\zkalm\Desktop\DendriteDetection\Data\',imNumber,'Blebs.tif')); 
