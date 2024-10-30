%%

clc;
% clearvars;
close all;

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

    noise(:,:,i) = img(:,:,i).*double(imcomplement(boundaries))./255;
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
set(gca,'FontSize',16)

figure
imshow(noise(:,:,1).') % Not sure if this is a proper method of tracking noise levels
title("Noisy Image from Last Week")
set(gca,'FontSize',16)

figure
imshow(img(:,:,1).')
title("Original Image")
set(gca,'FontSize',16)

%% Isolate noisy land pixels

denoise_land = zeros(size(img));
noisy_land = zeros(size(img));
boundaries_logical = imcomplement(boundaries)./255; % land = 1, water = 0

for i = 1:N
    disp(strcat(num2str(i), "/", num2str(N)))
    orig_denoise = wdenoise2(img(:,:,i));

    fig1 = figure;
    axis on;
    imshow(orig_denoise.');
    title(strcat("Interferogram #", num2str(i), " After wdenoise2()"))
    set(gca,'FontSize',16)

    single_denoise_land = orig_denoise.*double(boundaries_logical);

    figure;
    imshow(single_denoise_land.');
    title(strcat("Interferogram #", num2str(i), " Land Only After wdenoise2()"))
    set(gca,'FontSize',16)
    exportgraphics(fig1,strcat('C:\Users\mmpho\OneDrive - Washington University in St. Louis\Year 4\Capstone\Denoised Land Images\denoised_',num2str(i),'.png'))

    single_noisy_land = zeros(size(single_denoise_land,1),size(single_denoise_land,2));
    single_noisy_land((single_denoise_land < mean(single_denoise_land,"all")-(std(single_denoise_land,1,"all")/2)) & (boundaries_logical == 1)) = 1;

    fig2 = figure;
    imshow((single_noisy_land.*255).');
    title("Isolated Noisy Land Pixels")
    set(gca,'FontSize',16)
    exportgraphics(fig2,strcat('C:\Users\mmpho\OneDrive - Washington University in St. Louis\Year 4\Capstone\Denoised Land Images\noisy_land_',num2str(i),'.png'))

    denoise_land(:,:,i) = single_denoise_land;
    noisy_land(:,:,i) = single_noisy_land;
    close all;
end

%% Original code (only for igram #1)

denoise_test = wdenoise2(img(:,:,1));

figure;
axis on;
imshow(denoise_test.');
title("Interferogram #1 After wdenoise2()")
set(gca,'FontSize',16)

denoise_land = denoise_test.*double(boundaries_logical);

figure;
imshow(denoise_land.');
title("Interferogram #1 Land Only After wdenoise2()")
set(gca,'FontSize',16)

% Test 2 - this worked better than the foor loop, way more efficient
% Trying logical indexing in the for loop to remove test_row and test_col

single_noisy_land = zeros(size(denoise_land,1),size(denoise_land,2));
single_noisy_land((denoise_land < mean(denoise_land,"all")-(std(denoise_land,1,"all")/2)) & (boundaries_logical == 1)) = 1;

figure;
imshow((single_noisy_land.*255).');
title("Isolated Noisy Land Pixels")
set(gca,'FontSize',16)


%% Check values in denoised img vs original img

% Only look for land pixels that aren't the same value
[diff_row, diff_col] = find((denoise_land ~= im1) & (denoise_land ~= 0));

figure
histogram(diff_row)
title("Rows")

figure;
histogram(diff_col)
title("Cols")

%% Scrap Code area
%
% This took 7 years to run and didn't look too great
% for i = 980:size(denoise_land,1)
%     for j = 1:size(denoise_land,2)
%         if boundaries_logical(i,j) == 1
%             if denoise_land(i,j) < mean(denoise_land,"all")-(std(denoise_land,1,"all")/2)
%                 noisy_land(i,j) = 1;
%             end
%         end
%     end
%     disp(strcat("i = ", num2str(i)))
% end
%
% figure;
% imshow((noisy_land.*255).');
% title("Noisy land pixels???")

% See what happens if we do another denoise lol
% Update: it did nothing lol so commenting out
% denoise_test2 = wdenoise2(denoise_land);
%
% figure;
% imshow(denoise_test2.');
% title("Wavelet Denoise x2")

% figure;
% imshow(L.');
% title("L")

% noisy_land = zeros(size(denoise_land,1),size(denoise_land,2));
% [test_row, test_col] = find((denoise_land < mean(denoise_land,"all")-(std(denoise_land,1,"all")/2)) & (boundaries_logical == 1));
%
% for i = 1:length(test_row)
%     noisy_land(test_row(i),test_col(i)) = 1;
%     % disp(strcat("i = ", num2str(i), "  test_row = ", num2str(test_row(i)), " test_col = ", num2str(test_col(i))))
% end
%
% figure;
% imshow((noisy_land.*255).');
% title("Noisy land pixels???")
% [nm_row, nm_col] = find(noisy_land ~= noisy_land2);

% no_x_tix = 5;
% no_y_tix = 5;
% [pix, labels] = igram_latlong(subdir,no_x_tix,no_y_tix);
% no_x_tix = no_x_tix + 1;
% no_y_tix = no_y_tix + 1;
%
% addpath('C:\Users\mmpho\radar-lab\sp24')
% subdir = 'sbas_24';
% addpath(strcat('C:\Users\mmpho\sent_test\',subdir))
