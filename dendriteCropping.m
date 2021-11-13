function [dend_crop] = dendriteCropping(img_number) 
% This crops images to just get dendrites 
    %% GET IMAGES
    close all
    %% Load images and make grayscale for processing
    parent = pwd;
    image_names = dir(fullfile(parent,'*.tif'));
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
    imshow(im_dend)
    
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
end

