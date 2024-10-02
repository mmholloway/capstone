%% Set Limits
% tic
img = unw_phase.*mask;
loopnumber = 5;
size_limit = 10;
minPeakValue = 0;
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
% graycompsave2 = zeros(length(imgsbound), length(imgsbound(:,1))); % holder variable  for blur
% noiseTest1 = zeros(length(imgsbound), length(imgsbound(:,1))); % holder variable for 3D plots
% withCellCountSave = zeros(length(imgsbound), length(imgsbound(:,1))); % holder variable for including cell counts
% grayscaleshow = zeros(length(imgsbound), length(imgsbound(:,1))); % debugging matrix used to see what the final grayscale image looks like
% bound = zeros(length(imgsbound), length(imgsbound(:,1))); % 

peaks = zeros(length(imgsbound), length(imgsbound(:,1))); % Matrix that is used to find the condensates
% peaksaveList = zeros(length(imgsbound), length(imgsbound(:,1))); % Matrix used to save the full pixels of the condensates
peakLocations = zeros(length(imgsbound)*2,3); % Matrix to store the x, y, and number of peaks for each plane
% peakLocationAll = zeros(length(imgsbound)*2, 3);
% peakLocationMatrixAll = zeros(length(imgsbound), length(imgsbound(:,1)));
loc = zeros(length(imgsbound)*3, 2);

% num_peaks_max = 0; % Maximum number of peaks in a plane
% graytestsave = zeros(length(imgsbound), length(imgsbound(:,1))); % Holder matrix for the graytest matrix (Only if blur)
% graycompsave = zeros(length(imgsbound), length(imgsbound(:,1))); % Holder matrix for the graycomp matrix

for i=1:length(img(1,1,:))
    disp(i + "/" + length(img(1,1,:)))
%% Grayscale
%     peakLocations = zeros(length(imgsbound),3); % The locations of the peaks on the current peak
        peaksave = zeros(length(imgsbound), length(imgsbound(:,1))); % 
        peaksave1 = zeros(length(imgsbound), length(imgsbound(:,1))); % A holder variable for peaksave

        imgs = img(:,:,i); % The specific image being analized 
%     imgsdoub = im2double(imgs)*255;
%     if blur
%         for replay = 1:5
%         graytest = rgb2gray(imgs);
%             for x = 2:length(imgs)-1
%                 for y = 2:length(imgs(:,1))-1
%                     noise_counter = 0;
%                     if graytest(x,y) >= graytestHigh
%                         graytest(x,y) = 255;
%                     elseif graytest(x,y) <= graytestLow
%                         graytest(x,y) = 0;
%                     end
%                     if graytest(x,y) == 255
%                         if graytest(x+1,y+1) == 0 || graytest(x,y) == 255
%                             noise_counter = noise_counter+1;
%                         end
%                         if graytest(x+1,y) == 0 || graytest(x,y) == 255
%                             noise_counter = noise_counter+1;
%                         end
%                         if graytest(x+1,y-1) == 0 || graytest(x,y) == 255
%                             noise_counter = noise_counter+1;
%                         end
%                         if graytest(x,y+1) == 0 || graytest(x,y) == 255
%                             noise_counter = noise_counter+1;
%                         end
%                         if graytest(x,y-1) == 0 || graytest(x,y) == 255
%                             noise_counter = noise_counter+1;
%                         end
%                         if graytest(x-1,y+1) == 0 || graytest(x,y) == 255
%                             noise_counter = noise_counter+1;
%                         end
%                         if graytest(x-1,y) == 0 || graytest(x,y) == 255
%                             noise_counter = noise_counter+1;
%                         end
%                         if graytest(x-1,y-1) == 0 || graytest(x,y) == 255
%                             noise_counter = noise_counter+1;
%                         end
%                         if noise_counter == 8
%                             graytestsave(x,y) = 0;
%                         else
%                             graytestsave(x,y) = graytest(x,y);
%                         end
%                     else
%                         graytestsave(x,y) = graytest(x,y);
%                     end
%                 end
%             end
% 
%             graytest = graytestsave;
%         end
%         for replay = 1:3
%         for x=2:length(graytest)-1
%             for y=2:length(graytest(:,1))-1
%                 graytestsave(x,y) = (graytest(x+1,y+1) + graytest(x+1,y) + graytest(x+1,y-1) +...
%                                      graytest(x,y+1)   + graytest(x,y)   + graytest(x,y-1) +...
%                                      graytest(x-1,y+1) + graytest(x-1,y) + graytest(x-1,y-1))/9;
%                 if graytestsave(x,y)<=graytestLow
%                     graytestsave(x,y)=0;
%                 end
%             end
%         end
%         graytest = graytestsave;
%         grayscalebound = uint8(graytest);
%     %     end
%     else
%         grayscalebound = rgb2gray(imgs);
%     end
%     grayscale = rgb2gray(imgs); % Turn the image into Grayscale
        grayscale = imgs;
        graycomp = imcomplement(grayscale); % flip the grayscale bright to dark
%     for x = 1:length(grayscale)
%         for y = 1:length(graycomp(:,1))
%             if grayscale(x,y)>grayscaleMin % If the pixel is brighter than the maximum grayscale complement value, remove the pixel from calculation
%                 graycomp(x,y)=255;
%             end
%             if grayscale(x,y)<=noiseMin % If the pixel is darker than the minimum grayscale complement value, remove the pixel from calculation
%                 graycomp(x,y)=255;
%             end
%             
%         end
%     end    
%     grayscalenew = zeros(length(grayscale),length(grayscale(:,1)));
%     graycompnew = zeros(length(imgsbound), length(imgsbound(:,1)));
%     grayscalesave = zeros(length(grayscale),length(grayscale(:,1)));
%     

    
%% Boundary
        [B, L] = bwboundaries(grayscale, 'holes');
        imgHold = label2rgb(L);
        
        boundaries = rgb2gray(imgHold); % Create a grayscale image of the boundaries
    
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
            for x = 2:length(grayscale)-1
                for y = 2:length(grayscale(1,:))-1
                    if (y>2 && y<length(imgsbound)-1 && x>2 && x<length(imgsbound(:,1))-1)
                        if (boundaries(x,y+1) == 255 && boundaries(x,y-1) == 255 || boundaries(x-1,y) == 255 && boundaries(x+1,y) == 255 ||...
                            boundaries(x,y+2) == 255 && boundaries(x,y-2) == 255 || boundaries(x-2,y) == 255 && boundaries(x+2,y) == 255 ||...
                            boundaries(x+1,y+1) == 255 && boundaries(x-1,y-1) == 255 || boundaries(x-1,y+1) == 255 && boundaries(x+1,y-1) == 255)
                                boundaries(x,y) = 255;
                        end
                    end
                end
            end
        end
    
        % If the pixel is outside the cell boundaries set the value to 255
        % (remove from calculations)
        for x=1:length(grayscale)
            for y=1:length(grayscale(1,:))
                if boundaries(x,y)==0
                    graycomp(x,y)=graycomp(x,y);
                else
                    graycomp(x,y)=255;
                end
            end
        end
      
%     if blur
%         for q=1:loopNumber
%             for x = 1:length(grayscale)
%                 for y = 1:length(grayscale(1,:))
%                     if bound(x,y) == 0
%                         if x >= 5 && x <= length(grayscale)-5
%                             if y >= 5 && y <= height(grayscale)-5
%                                 if graycomp(x,y)>0
%                                     graycompsave2(x,y) = ((graycomp(x+1,y+1) + graycomp(x+1,y) + graycomp(x+1,y-1) + ...
%                                                          graycomp(x,y+1) + graycomp(x,y-1) + ...
%                                                          graycomp(x-1,y+1) + graycomp(x-1,y) + graycomp(x-1,y-1))/8+graycomp(x,y))/2;
%                                 else
%                                     graycompsave2(x,y) = (graycomp(x+1,y+1) + graycomp(x+1,y) + graycomp(x+1,y-1) + ...
%                                                           graycomp(x,y+1) + graycomp(x,y-1) + graycomp(x,y) +...
%                                                           graycomp(x-1,y+1) + graycomp(x-1,y) + graycomp(x-1,y-1))/9;
%                                 end
%                             end
%                         end
%                     else
%                         graycompsave2(x,y) = 255;
%                     end
%                 end
%             end
%             graycomp = graycompsave2;
%         end
%         for q=1:2
%             for x=2:length(graycomp)-1
%                 for y = 2:length(graycomp(:,1))-1
%                     blackcount = 0;
%                     if graycomp(x+1,y) == 0
%                         blackcount = blackcount+1;
%                     end
%                     if graycomp(x,y+1) == 0
%                         blackcount = blackcount+1;
%                     end
%                     if graycomp(x-1,y) == 0
%                         blackcount = blackcount+1;
%                     end
%                     if graycomp(x,y-1) == 0
%                         blackcount = blackcount+1;
%                     end
%                     if blackcount<2
%                         graycomp(x,y) = 255;
%                     end
%                 end
%             end
%         end
%     end

%% Call Find_Peaks
%     tic
%     [peaks(:,:), loc(:,:), Num_Peaks, peakLocMatrix] = Find_Peaks(imcomplement(graycomp), 10, size_limit, minPeakValue);
%     toc
%     % create a holder location list that also holds the number of peaks
%     for x = 1:Num_Peaks 
%         for y = 1:2 
%             peakLocations(x,y)=loc(x,y);
%         end
%         peakLocations(x,3) = Num_Peaks;
%     end
%     for list=1:Num_Peaks
%         peaksave(loc(list,1),loc(list,2)) = 255;
%     end
%     
%     for q=1:5
%         for x = 3:length(peaks)-3
%             for y = 3:length(peaks(:,1))-3
%                 if peaksave(x,y) == 255 || peaksave(x+1,y) == 255 && graycomp(x,y) <=grayscaleMin || peaksave(x-1,y) == 255 && graycomp(x,y) <=grayscaleMin ||...
%                    peaksave(x,y+1) == 255 && graycomp(x,y) <=grayscaleMin || peaksave(x,y-1) == 255 && graycomp(x,y) <=grayscaleMin || peaksave(x+2,y) == 255 && graycomp(x,y) <=grayscaleMin ||...
%                    peaksave(x-2,y) == 255 && graycomp(x,y) <=grayscaleMin || peaksave(x,y+2) == 255 && graycomp(x,y) <=grayscaleMin || peaksave(x,y-2) == 255 && graycomp(x,y) <=grayscaleMin
%                     peaksave1(x,y) = 255;
%                 end
%             end
%         end
%         peaksave = peaksave1;
%     end
    
%% Expand Out
%     while change==1 % while there is a pixel added to the size of a condensate
%         changeflag = 0; % Flag to show if there has been a change in condensate size
%         q=q+1; % Counter only for timing
%         disp(q)
% %         % Expand out the peak to fill the area as shown in Graycompsave, this
% %         % is done by adding the two nearest pixels in the four cardinal
% %         % directions to each condensate and repeating the process until no
% %         % pixels are added
%         for x = 3:length(peaks)-3
%             for y = 3:length(peaks(:,1))-3
%                 for z = 1:lastIMG
%                     if peaksave(x,y)>0
%                         withCellCountSave(x,y) = peaksave(x,y);
%                     elseif peaksave(x+1,y)>0 && graycomp(x,y)<255
%                         withCellCountSave(x,y) = peaksave(x+1,y);
%                         changeflag = 1; % A flag showing that a pixel was added
%                         changecount = changecount+1;
%                     elseif peaksave(x+2,y)>0 && graycomp(x,y)<255  
%                         withCellCountSave(x,y) = peaksave(x+2,y);
%                         changeflag = 1;
%                         changecount = changecount+1;
%                     elseif peaksave(x,y+1)>0 && graycomp(x,y)<255  
%                         withCellCountSave(x,y) = peaksave(x,y+1);
%                         changeflag = 1;
%                         changecount = changecount+1;
%                     elseif peaksave(x,y+2)>0 && graycomp(x,y)<255 
%                         withCellCountSave(x,y) = peaksave(x,y+2);
%                         changeflag = 1;
%                         changecount = changecount+1;
%                     elseif peaksave(x-1,y)>0 && graycomp(x,y)<255   
%                         withCellCountSave(x,y) = peaksave(x-1,y);
%                         changeflag = 1;
%                         changecount = changecount+1;
%                     elseif peaksave(x-2,y)>0 && graycomp(x,y)<255   
%                         withCellCountSave(x,y) = peaksave(x-2,y);
%                         changeflag = 1;
%                         changecount = changecount+1;
%                     elseif peaksave(x,y-1)>0 && graycomp(x,y)<255   
%                         withCellCountSave(x,y) = peaksave(x,y-1);
%                         changeflag = 1;
%                         changecount = changecount+1;
%                     elseif peaksave(x,y-2)>0 && graycomp(x,y)<255   
%                         withCellCountSave(x,y) = peaksave(x,y-2);
%                         changeflag = 1;
%                         changecount = changecount+1;
%                     end
%                 end
%             end
%         end
%         peaksave = withCellCountSave;
%         if ~changeflag
%             change = 0;
%         else
%             change = 1;
%         end
%     end
% %% Full Matricies
%     peaksaveAll(:,:,i) = peaksave(:,:);
%     graycompAll(:,:,i) = graycomp(:,:);
%     boundariesAll(:,:,i) = boundaries(:,:);
end
%% Figures
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

%% Functions
% function array_changed = filter_through_x(img,k)
%     parfor i = 1:length(img)
%         img_changed(i,:,k) = filter_through_y(img,k,i)
%     end
%     disp(size(img_changed))
% end
% 
% function array_changed = filter_through_y(img,k,i)
%     parfor j = 1:length(img(1,:,1))
%         if isnan(img(i,j,k))
%             img(i,j,k) = 0;
%         end
%     end
%     array_changed = img(i,:,k);
% end