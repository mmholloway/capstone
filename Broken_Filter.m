close all
%% Set Limits
% tic
img = unw_phase.*mask;
loopnumber = 5;
plus_minus = 2;
% length(img(1,1,:))
% tic
% parfor k = 1:51
%     img_changed(:,:,k) = filter_through_x(img,k);
% end
% toc
for k = 1:length(img(1,1,:))
    for j = 1:length(img(1,:,1))
        for i = 1:length(img(:,1,1))
            if isnan(img(i,j,k))
                img(i,j,k) = 0;
            end
        end
    end
end
% for j = 1:length(img(1,:,1))
%     for i = 1:length(img(:,1,1))
%         img_average(i,j)
%     end
% end

% img = img_changed;
% figure
% imshow(transpose(img(:,:,1)))




%% Initilize variable
imgsbound = img; % Used to determine the number of pixels in the image


peaks = zeros(length(imgsbound), length(imgsbound(:,1))); % Matrix that is used to find the condensates

peakLocations = zeros(length(imgsbound)*2,3); % Matrix to store the x, y, and number of peaks for each plane

loc = zeros(length(imgsbound)*3, 2);

for i=1:length(img(1,1,:))
    clear pixel_values
% for i = 1:2
    disp(i + "/" + length(img(1,1,:)))
%% Grayscale
%     peakLocations = zeros(length(imgsbound),3); % The locations of the peaks on the current peak
        peaksave = zeros(length(imgsbound), length(imgsbound(:,1))); % 
        peaksave1 = zeros(length(imgsbound), length(imgsbound(1,:))); % A holder variable for peaksave
        boundaries_save = ones(length(imgsbound),length(imgsbound(1,:,1)))*255;
        imgs = img(:,:,i); % The specific image being analized 
        grayscale = imgs;
        graycomp = imcomplement(grayscale); % flip the grayscale bright to dark
    
%% Boundary
        [B, L] = bwboundaries(grayscale, 'holes');
        imgHold = label2rgb(L);
        
        boundaries = rgb2gray(imgHold); % Create a grayscale image of the boundaries
        
        for j = 1:250
            pixel_values(j) = sum(sum(boundaries == j));
        end
%         figure
%         plot(pixel_values)
        pixel_lim = find(pixel_values == max(pixel_values));
        pixel_lim_min = pixel_lim-plus_minus;
        pixel_lim_max = pixel_lim+plus_minus;

        % The First and Last rows and columns removed from the boundaries and
        % calculations (would be removed later)
        boundaries(:,1) = 255;
        boundaries(:,2) = 255;
        boundaries(:,length(grayscale(1,:))-1) = 255;
        boundaries(:,length(grayscale(1,:))) = 255;
        boundaries(1,:) = 255;
        boundaries(2,:) = 255;
        boundaries(height(grayscale)-1, :) = 255;
        boundaries(height(grayscale), :) = 255;
    
        % If the boundary pixel is brighter than the limit, remove it from the
        % boundary, otherwhise set it to 0
%         for x = 3:length(grayscale)-2
%             for y = 3:length(grayscale(1,:))-2
%     %           Finalize Boundary of Cell
%                 if (boundaries(x,y)>boundaryLimit) % how dark of a cell wall is removed
%                     boundaries(x,y) = 255;
%                 else
%                     boundaries(x,y) = 0;
%                 end
%             end
%         end
        
        % noise removal from boundaries
        for loop = 1:loopnumber
            for x = 3:length(grayscale)-2
                for y = 3:length(grayscale(1,:))-2
                    if (y>2 && y<length(imgsbound)-1 && x>2 && x<length(imgsbound(:,1))-1)
%                         if boundaries(x,y)>pixel_lim_max || boundaries(x,y)<pixel_lim_min
%                             boundaries_save(x,y)=255;
%                         else
%                             boundaries_save(x,y) = pixel_lim;
%     %                         disp("hi")
%                         end
%                         boundaries = boundaries_save;

                        if (boundaries(x,y+1)<pixel_lim_min && boundaries(x,y-1)<pixel_lim_min || boundaries(x-1,y)<pixel_lim_min && boundaries(x+1,y)<pixel_lim_min ||...
                                boundaries(x,y+2)<pixel_lim_min && boundaries(x,y-2)<pixel_lim_min || boundaries(x-2,y)<pixel_lim_min && boundaries(x+2,y)<pixel_lim_min ||...
                                boundaries(x+1,y+1)<pixel_lim_min && boundaries(x-1,y-1)<pixel_lim_min || boundaries(x-1,y+1)<pixel_lim_min && boundaries(x+1,y-1)<pixel_lim_min)
                                    boundaries_save(x,y) = 255;
                        
                        elseif (boundaries(x,y+1)<=pixel_lim_max && boundaries(x,y-1)<=pixel_lim_max || boundaries(x-1,y)<=pixel_lim_max && boundaries(x+1,y)<=pixel_lim_max ||...
                                boundaries(x,y+2)<=pixel_lim_max && boundaries(x,y-2)<=pixel_lim_max || boundaries(x-2,y)<=pixel_lim_max && boundaries(x+2,y)<=pixel_lim_max ||...
                                boundaries(x+1,y+1)<=pixel_lim_max && boundaries(x-1,y-1)<=pixel_lim_max || boundaries(x-1,y+1)<=pixel_lim_max && boundaries(x+1,y-1)<=pixel_lim_max)
                                    boundaries_save(x,y) = pixel_lim;
%                                     disp(x + "," + y)
                        else
                            boundaries_save(x,y) = 255;
                        end
                    end
                end
            end
            
            boundaries = uint8(boundaries_save);
        end
    
        % If the pixel is outside the cell boundaries set the value to 255
        % (remove from calculations)
        for x=1:length(grayscale)
            for y=1:length(grayscale(1,:))
                if boundaries(x,y)==pixel_lim
                    graycomp(x,y)=1;
                else
                    graycomp(x,y)=0;
                end
            end
        end

noise(:,:,i) = abs(img(:,:,i)-double(boundaries)).*graycomp;
    for x = 1:length(noise)
        for y = 1:length(noise(1,:,1))
            if noise(x,y,i) == pixel_lim
                noise(x,y,i) = 1;
            else
                noise(x,y,i) = 0;
            end
        end
    end
end
%% Figures
figure
imshow(transpose(boundaries))
title("Filtered")
figure
imshow(transpose(noise(:,:,2)))
title("noise")
figure
imshow(transpose(img(:,:,2)))
title("img")
%     figure
%     imshow(uint8(grayscale))
%     title("Yay")
% figure
% imshowpair(peaksave(:,:,51), graycomp,'montage')
% title("Final Check")
% toc
%     imshowpair(peaks(:,:,i),graycomp)
%     figure
%     imshow(boundaries)
%     title("Boundaries")
%     figure
%     imshow(imgs)
