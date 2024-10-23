%% Set Limits
% tic
img = unw_phase.*mask;
loopnumber = 3;
size_limit = 10;
minPeakValue = 0;
pixel_lim = 40;
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

%% Initilize variable
imgsbound = img; % Used to determine the number of pixels in the image

% peaks = zeros(length(imgsbound), length(imgsbound(:,1))); % Matrix that is used to find the condensates
% peakLocations = zeros(length(imgsbound)*2,3); % Matrix to store the x, y, and number of peaks for each plane
% loc = zeros(length(imgsbound)*3, 2);

boundaries_mtx = zeros(size(img));

for i=1:length(img(1,1,:))
    disp(i + "/" + length(img(1,1,:)))
    % Grayscale
    peaksave = zeros(length(imgsbound), length(imgsbound(:,1))); %
    peaksave1 = zeros(length(imgsbound), length(imgsbound(:,1))); % A holder variable for peaksave

    imgs = img(:,:,i); % The specific image being analized
    grayscale = imgs;
    graycomp = imcomplement(grayscale); % flip the grayscale bright to dark


    % Boundary
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

    % noise removal from boundaries
    for loop = 1:loopnumber
        for x = 3:length(grayscale)-2
            for y = 3:length(grayscale(1,:))-2
                if (y>2 && y<length(imgsbound)-1 && x>2 && x<length(imgsbound(:,1))-1)
                    if (boundaries(x,y+1)>pixel_lim && boundaries(x,y-1)>pixel_lim || boundaries(x-1,y)>pixel_lim && boundaries(x+1,y)>pixel_lim ||...
                            boundaries(x,y+2)>pixel_lim && boundaries(x,y-2)>pixel_lim || boundaries(x-2,y)>pixel_lim && boundaries(x+2,y)>pixel_lim ||...
                            boundaries(x+1,y+1)>pixel_lim && boundaries(x-1,y-1)>pixel_lim || boundaries(x-1,y+1)>pixel_lim && boundaries(x+1,y-1)>pixel_lim)
                        boundaries(x,y) = 255;
                    end
                    if (boundaries(x,y+1)<pixel_lim && boundaries(x,y-1)<pixel_lim || boundaries(x-1,y)<pixel_lim && boundaries(x+1,y)<pixel_lim ||...
                            boundaries(x,y+2)<pixel_lim && boundaries(x,y-2)<pixel_lim || boundaries(x-2,y)<pixel_lim && boundaries(x+2,y)<pixel_lim ||...
                            boundaries(x+1,y+1)<pixel_lim && boundaries(x-1,y-1)<pixel_lim || boundaries(x-1,y+1)<pixel_lim && boundaries(x+1,y-1)<pixel_lim)
                        boundaries(x,y) = 0;
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

    noise(:,:,i) = img(:,:,i)-double(imcomplement(boundaries));
    boundaries_mtx(:,:,i) = boundaries;

end
%% Figures

% for i = 1:size(boundaries_mtx,3)
%     figure
%     imshow(transpose(boundaries_mtx(:,:,i)))
%     title(strcat("Interferogram ", num2str(i)));
%     exportgraphics(gcf,strcat('C:\Users\mmpho\OneDrive - Washington University in St. Louis\Year 4\Capstone\IG Error Plots\Boundaries no', num2str(i),'.png'))
%     close;
% end

% for i = 1:size(boundaries_mtx,3)
%     for j = 1:size(boundaries_mtx,3)
%         if ~(sum(boundaries_mtx(:,:,i) == boundaries_mtx(:,:,j),"all") == 1890*1016)
%             disp(strcat("Not matching - i = ", num2str(i), " / j = ", num2str(j)))
%         end
%     end
% end

figure
imshow(transpose(boundaries))
title("Boundaries");
figure
imshow(transpose(noise(:,:,2)))
title("noise")
figure
imshow(transpose(img(:,:,2)))
title("img")

%% Isolate noisy land pixels

