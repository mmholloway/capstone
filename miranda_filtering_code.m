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
den_coeffs = [];
orig_coeffs = [];
ss = [];
shifts = [];

for i = 1:N
    disp(strcat(num2str(i), "/", num2str(N)))
    [orig_denoise,den_coeff,orig_coeff,s,shift] = wdenoise2(img(:,:,i));
    % den_coeffs = [den_coeffs; den_coeff];
    % orig_coeffs = [orig_coeffs; orig_coeff];
    % ss = [ss; s];
    % shifts = [shifts; shift];

    % fig1 = figure;
    % axis on;
    % imshow(orig_denoise.');
    % title(strcat("Interferogram #", num2str(i), " After wdenoise2()"))
    % set(gca,'FontSize',16)

    single_denoise_land = orig_denoise.*double(boundaries_logical);

    % fig1 = figure;
    % imshow(single_denoise_land.');
    % title(strcat("Interferogram #", num2str(i), " Land Only After wdenoise2()"))
    % set(gca,'FontSize',16)
    % exportgraphics(fig1,strcat('C:\Users\mmpho\OneDrive - Washington University in St. Louis\Year 4\Capstone\Denoised Land Images\BRISQUE Trainers\brisque_model',num2str(i),'.png'))

    single_noisy_land = zeros(size(single_denoise_land,1),size(single_denoise_land,2));
    single_noisy_land((single_denoise_land < mean(single_denoise_land,"all")-(std(single_denoise_land,1,"all")/2)) & (boundaries_logical == 1)) = 1;

    % fig2 = figure;
    % imshow((single_noisy_land.*255).');
    % title("Isolated Noisy Land Pixels")
    % set(gca,'FontSize',16)
    % exportgraphics(fig2,strcat('C:\Users\mmpho\OneDrive - Washington University in St. Louis\Year 4\Capstone\Denoised Land Images\noisy_land_',num2str(i),'.png'))

    denoise_land(:,:,i) = single_denoise_land;
    noisy_land(:,:,i) = single_noisy_land;
    % close all;
end

%% Original code (only for igram #1)

[denoise_test,igram1_denoised_coeffs,igram1_org_coeffs,s,shifts] = wdenoise2(img(:,:,1));

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

figure;
tiledlayout(2,2)

nexttile
imshow(img(:,:,1).')
title("Interferogram #1")
set(gca,'FontSize',16)

nexttile
imshow(denoise_test.');
title("Interferogram #1 After wdenoise2()")
set(gca,'FontSize',16)

nexttile
imshow((single_noisy_land.*255).');
title("Isolated Noisy Land Pixels")
set(gca,'FontSize',16)

nexttile
imshow(denoise_land.');
title("Interferogram #1 Land Only After wdenoise2()")
set(gca,'FontSize',16)

%% Looking at wavelet coefficients

% same_coeffs = igram1_denoised_coeffs == igram1_org_coeffs;
% disp(sum(same_coeffs))

% same_idx = find(same_coeffs);
% idx_jumps = diff(same_idx);
% jumps_loc = find(idx_jumps > 1);
% jumps_loc = [jumps_loc; find(idx_jumps > 10), zeros(1,length(jumps_loc) - length(find(idx_jumps > 10)))];

% figure;
% hold on
% histogram(jumps_loc(1,:))
% histogram(jumps_loc(2,:))
% legend("> 1", "> 10")
%
% figure;
% plot(igram1_denoised_coeffs)

% figure;
% noise_pixels = single_noisy_land.*img(:,:,1);
% flat_snl = reshape(noise_pixels,1,[]);
% plot(fft(flat_snl))

good_den_coeffs = den_coeffs(19,:);
bad_den_coeffs = den_coeffs(49,:);

good_orig_coeffs = orig_coeffs(19,:);
bad_den_coeffs = orig_coeffs(49,:);

figure;
tiledlayout(1,2)

nexttile
hold on;
plot(good_den_coeffs)
plot(bad_den_coeffs)
title("Denoised Coefficient Comparison");
legend("IG #19", "IG #49")
xlabel("Coefficient Number")
ylabel("Coefficient Value")
set(gca,'FontSize',16)
axis tight

nexttile
hold on;
plot(good_den_coeffs)
plot(bad_den_coeffs)
title("Original Coefficient Comparison");
legend("IG #19", "IG #49")
xlabel("Coefficient Number")
ylabel("Coefficient Value")
set(gca,'FontSize',16)
axis tight

%% Calculate Contrast-to-Noise ratio (CNR)
% https://en.wikipedia.org/wiki/Contrast-to-noise_ratio
% Calculate per the information here: https://link.springer.com/content/pdf/10.1007/978-0-387-73507-8_7.pdf
% Let signal intensity be the average of the denoised land image, masked by
% the boundaries matrix
% Let the noise intensity be the average of the isolated noise pixels
% Then, the CNR is these two averages divided by the standard deviation of
% the image as calculated in the second link

cnr = zeros(1,size(unw_phase,3));

for i = 1:N
    s_land = mean(nonzeros(denoise_land(:,:,i).*double(boundaries)),"all");
    s_noise = mean(nonzeros(noisy_land(:,:,i).*unw_phase(:,:,i)),"all");

    stdev_land = std(denoise_land(:,:,i).*double(boundaries),0,"all");
    stdev_noise = std(noisy_land(:,:,i).*unw_phase(:,:,i),0,"all");
    stdev = (1/2)*sqrt(var(denoise_land(:,:,i).*double(boundaries),0,"all")+var(noisy_land(:,:,i).*unw_phase(:,:,i),0,"all"));

    cnr(i) = (s_land-s_noise)/stdev_noise;
end

good_resp = [1, 5, 7, 15:17, 18, 19, 21, 23, 24, 27:30, 32:34, 41, 43, ...
    44, 47, 48, 50, 51];
bad_resp = [2:4, 6, 8:14, 20, 22, 25, 26, 31, 35:40, 42, 45, 46, 49];

for i = 1:length(good_resp)
    disp(strcat("GOOD interferogram idx = ", num2str(good_resp(i)), ...
        " CNR = ", num2str(cnr(good_resp(i)))))
end

for i = 1:length(bad_resp)
    disp(strcat("BAD interferogram idx = ", num2str(bad_resp(i)), ...
        " CNR = ", num2str(cnr(bad_resp(i)))))
end

figure;
tiledlayout(1,2)

nexttile
boxplot(cnr(good_resp))
title("CNR of Interferograms Responding Well to wdenoise2()")
xlabel("Good Interferograms")
ylabel("CNR Value")
set(gca,'FontSize',16)

nexttile
boxplot(cnr(bad_resp))
title("CNR of Interferograms Responding Poorly to wdenoise2()")
xlabel("Bad Interferograms")
ylabel("CNR Value")
set(gca,'FontSize',16)

%% Evaluate more image quality metrics

% Implement BRISQUE
brisque_scores = zeros(1,size(unw_phase,3));
for i = 1:length(brisque_scores)
    denoise_mask = denoise_land(:,:,i) ~= 0;
    img = unw_phase(:,:,i).*denoise_mask;
    min_val = min(img(:));
    max_val = max(img(:));
    img_scaled = (img - min_val) / (max_val - min_val);

    brisque_scores(i) = brisque(img_scaled);
end

figure;
boxplot(brisque_scores)
title("BRISQUE Scores", "No Model")
xlabel("Interferograms")
ylabel("BRISQUE Score")

figure;
plot(brisque_scores)
title("BRISQUE Scores", "No Model")
xlabel("Interferograms")
ylabel("BRISQUE Score")

%% Implement BRISQUE with a model of 5 images
imds = imageDatastore("C:\Users\mmpho\OneDrive - Washington University in St. Louis\Year 4\Capstone\Denoised Land Images\BRISQUE Trainers");
% opinion_scores = [90, 55, 95, 10, 10]; % produced some bad results
% opinion_scores = [10, 30, 2, 100, 100];
% opinion_scores = [0 50 0 100 100]/4; % with Kam input
opinion_scores = [0 50 0 30 100 100 100 0 100 0]/2; % with Kam input
brisque_model = fitbrisque(imds,opinion_scores);

bs_model = zeros(1,size(unw_phase,3));
for i = 1:N
    bs_model(i) = brisque(denoise_land(:,:,i), brisque_model);
end

figure;
boxplot(bs_model)
title("BRISQUE Scores", "With Model")
xlabel("Interferograms")
ylabel("BRISQUE Score")

figure;
plot(bs_model)
title("BRISQUE Scores", "With Model")
xlabel("Interferograms")
ylabel("BRISQUE Score")

%% Implement PIQE model
piqe_scores = zeros(1,size(unw_phase,3));
for i = 1:length(piqe_scores)
    piqe_scores(i) = piqe(denoise_land(:,:,i));
end

figure;
boxplot(piqe_scores)
title("PIQE Scores")
xlabel("Interferograms")
ylabel("PIQE Score")

figure;
plot(piqe_scores)
title("PIQE Scores")
xlabel("Interferograms")
ylabel("PIQE Score")

%% Plot images for BRISQUE model

for i = [3 10 13 41 44]
    fig = figure;
    imshow((denoise_land(:,:,i)).');
    title(strcat("Interferogram #", num2str(i)));
    set(gca,'FontSize',16)
    exportgraphics(fig,strcat('C:\Users\mmpho\OneDrive - Washington University in St. Louis\Year 4\Capstone\Denoised Land Images\BRISQUE Trainers\brisque_model',num2str(i),'.png'))
    close all;
end

%% Plotting wrapped phase per Roger's recommendation

% Import wrapped phase files
subdir = 'sbas_24';
addpath(strcat('C:\Users\mmpho\sent_test\',subdir))

nr=1890; naz=1016;
wr_phase=zeros(nr,naz,N);

for i = 1:N
    disp(i + "/" + N)
    % wrapped phase
    filename=sprintf('%s',cells{i});
    fid=fopen(filename);
    dat=fread(fid,[2*nr,inf],'float','ieee-le');
    temp=dat(1:2:2*nr-1,:)+1i*dat(2:2:2*nr,:);
    wr_phase(:,:,i)=temp;
    fclose(fid);
end

%% Plot histogram of "noisy" and "non-noisy" land pixel phase values

wrapped_noise = wr_phase.*noisy_land;
denoise_land_log = double(denoise_land ~= 0);
wrapped_land = wr_phase.*denoise_land_log;

for i = 14
    % "Noisy" phase values
    fig1 = figure;
    tiledlayout(2,1)
    set(gcf, 'Position', get(0, 'Screensize'));

    nexttile
    histogram(angle(nonzeros(wrapped_noise(:,:,i))))
    axis tight
    title("Noisy Land Pixel Wrapped Phase Distribution", strcat("Interferogram #",num2str(i)))
    xlabel("Wrapped Phase Values (radians)")
    ylabel("Number of Pixels")
    % set(gca,'FontSize',16)
    % exportgraphics(fig1,strcat('C:\Users\mmpho\OneDrive - Washington University in St. Louis\Year 4\Capstone\Wrapped Phase Distributions\noisy_pix_',num2str(i),'.png'))

    % "Good" phase values
    nexttile
    % fig2 = figure;
    % set(gcf, 'Position', get(0, 'Screensize'));
    histogram(angle(nonzeros(wrapped_land(:,:,i))))
    axis tight
    title("Non-noisy Land Pixel Wrapped Phase Distribution", strcat("Interferogram #",num2str(i)))
    xlabel("Wrapped Phase Values (radians)")
    ylabel("Number of Pixels")
    % set(gca,'FontSize',16)
    exportgraphics(fig1,strcat('C:\Users\mmpho\OneDrive - Washington University in St. Louis\Year 4\Capstone\Wrapped Phase Distributions\comparison',num2str(i),'.png'))
    close all;
end

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
%
% net = denoisingNetwork("DnCNN");
% denoisedI = denoiseImage(img(:,:,1),net);
%
% figure;
% imshow(denoisedI);
%
%[diff_row, diff_col] = find((denoise_land ~= im1) & (denoise_land ~= 0));
%
% figure
% histogram(diff_row)
% title("Rows")
%
% figure;
% histogram(diff_col)
% title("Cols")


