%
% EPSc 407 Spring 2023
%
% This mfile provides skeleton code to perform  image classification
% of a 7-band MODIS scene.
%
clearvars, close all, clc;
%
%% PROBLEM 1) READ IN THE IMAGE AND LOOK AT THE BANDS
% read in and look at all 7 bands with imagesc
%

red  = uint8(load('modis1.dat')); % red
nir  = uint8(load('modis2.dat')); % near infrared
blue  = uint8(load('modis3.dat')); % blue
mir1 = uint8(load('modis4.dat')); % green
mir2 = uint8(load('modis5.dat')); % mid infrared
mir3  = uint8(load('modis6.dat'));
mir4 = uint8(load('modis7.dat'));

%
% look at each band separately and an rgb image
%

figure;
tiledlayout(2,4)

nexttile, imshow(cat(3,red,mir1,blue)), title("Composite Image");
nexttile, imshow(red), colorbar('horiz'), title("Red Image");
nexttile, imshow(nir), colorbar('horiz'), title("NIR Image");
nexttile, imshow(blue), colorbar('horiz'), title("Blue Image");
nexttile, imshow(mir1), colorbar('horiz'), title("Green Image");
nexttile, imshow(mir2), colorbar('horiz'), title("MIR Band 5 Image");
nexttile, imshow(mir3), colorbar('horiz'), title("MIR Band 6 Image");
nexttile, imshow(mir4), colorbar('horiz'), title("MIR Band 7 Image");

%% CREATE A TRAINING SET
% Identify the location (row,column) of "training pixels" and assign them a "group number":
% You will need to identify at least 2 training pixels per group
%

% Each line of tpix should have 3 numbers in the following format:
%     row,col,group
%      ^   ^  ^

% Group 1 - melt ponds
% Group 2 - clouds
% Group 3 - Charcot Island
% Group 4 - region of disintegration
% Group 5 - sea ice
% Group 6 - Lataday Island
% Group 7 - Wilkins Ice Shelf
% Group 8 - outside image

% tpix = [283, 42, 1;
%     288, 17, 1;
%     280, 39, 1;
%     376, 277, 2;
%     305, 212, 2;
%     97, 93, 2;
%     131, 157, 2;
%     120, 188, 3;
%     67, 190, 3;
%     182, 186, 4;
%     210, 234, 4;
%     182, 220, 4;
%     122, 272, 5;
%     155, 298, 5;
%     131, 320, 5;
%     159, 389, 6;
%     182, 324, 6;
%     138, 357, 6;
%     255, 202, 7;
%     284, 75, 7;
%     383, 137, 7;
%     253, 209, 7;
%     239, 115, 7];

tpix = [289, 19, 1;
        282, 41, 1;
        332, 270, 2;
        123, 129, 2;
        101, 186, 3;
        99, 155, 3;
        164, 193, 4;
        200, 210, 4;
        134, 309, 5;
        158, 319, 5;
        183, 358, 6;
        213, 278, 6;
        282, 96, 7;
        221, 169, 7;
        64, 28, 8;
        41, 269, 8;
        194, 41, 8];

row=tpix(:,1);   % y-value
col=tpix(:,2);   % x-value
group=tpix(:,3); % group number
ngroup=max(group);

%
% From these pixels, make a training set consisting of each training pixel's band values
% "train" should have the same number of rows as the number of training pixels, and the
% same number of columns as number of bands (in the MODIS case, 7).
%

train = zeros(length(tpix),7);
train_idx = 1;
for i = 1:size(row,1)
    train(train_idx,1:7) = [red(row(i), col(i)), nir(row(i), col(i)), blue(row(i), col(i)),...
        mir1(row(i), col(i)), mir2(row(i), col(i)), mir3(row(i), col(i)),...
        mir4(row(i), col(i))];
    train_idx = train_idx + 1;
end

%
%% PROBLEM 3A) PREPARE DATA FOR CLASSIFICATION
% Reshape image into one long vector of pixel band values.
% Convert from uint8 to double for classification.
%
nx=400;
ny=400;
N=nx*ny;
nz=7;

all_pix = [reshape(red,[N,1]), reshape(nir,[N,1]), reshape(blue,[N,1]),...
    reshape(mir1,[N,1]), reshape(mir2,[N,1]), reshape(mir3,[N,1]),...
    reshape(mir4,[N,1])];

all_pix = double(all_pix);
train = double(train);
%
%% PROBLEM 3B) CLASSIFY THE IMAGE
% Classify the image. Matlab's "classify" function requires the Statistics toolbox.
%   train and sample must have same number of columns.
%   train and group must have same number of rows.
%   misfit is nx-by-ny-by-ngroups and has probability (0-1) that each is a member of that group
% This may take up to a minute.  tic and toc will time the calculation
%

tic
%
% Perform classification
%
[class,err,misfit] = classify(all_pix,train,group);

toc

%
% PROBLEM 3C) LOOK AT THE CLASSIFIED IMAGE
% Reshape the Class vector and Each group error vetor (Nx1) back into an nx x ny matrix.
%

class = reshape(class,nx,ny);
misfit = reshape(misfit,nx,ny,ngroup);

%
% Visualize the classification.
% You could do this with any colormap, or you can make your own with RGB values
% for each Group Number that make sense to you. Feel free to use the map below
% or make your own.
%
% map=[1,0.96863,0.92157;...          % Group 1: Clouds
%     0.72941,0.83137,0.95686;...    % Group 2: Sea Ice
%     0,0.74902,0.74902;...          % Group 3: Floating Ice Shelf
%     0.95294,0.87059,0.73333;       % Group 4: Land Ice
%     0.078431,0.16863,0.54902];     % Group 5: Open Ocean

map = 1/255.*[39, 64, 1; 130, 138, 0; 242, 159, 5; 242, 92, 5; 214, 86, 140; 77, 133, 132; 166, 47, 3; 64, 13, 1]

figure;
colormap(map), colorbar;
image(class);
title("Classified Image","1 Melt Ponds, 2 Clouds, 3 Charcot Island, 4 Disintegration, 5 Sea Ice, 6 Lataday Island, 7 Wilkins Ice Shelf");

%
%% PROBLEM 3) OPTIONAL: LOOK AT THE MISFIT
%

figure;
tiledlayout(2,4);
nexttile, imagesc(misfit(:,:,1)), colorbar, title("Melt Ponds");
nexttile, imagesc(misfit(:,:,2)), colorbar, title("Clouds");
nexttile, imagesc(misfit(:,:,3)), colorbar, title("Charcot Island");
nexttile, imagesc(misfit(:,:,4)), colorbar, title("Disintegration");
nexttile, imagesc(misfit(:,:,5)), colorbar, title("Sea Ice");
nexttile, imagesc(misfit(:,:,6)), colorbar, title("Lataday Island");
nexttile, imagesc(misfit(:,:,7)), colorbar, title("Wilkins Ice Shelf");

M=max(misfit,[],3);
figure;
imagesc(M),colorbar,title('maximum fit (1=perfect,0=not at all)');

I = find(M < 0.90);
class2 = class;
class2(I) = ngroup + 1;

map = 1/255.*[39, 64, 1; 130, 138, 0; 242, 159, 5; 242, 92, 5; 214, 86, 140; 77, 133, 132; 166, 47, 3; 64, 13, 1; 0, 0, 0];

figure;
colormap(map), colorbar;
image(class);
title("Classified Image","1 Melt Ponds, 2 Clouds, 3 Charcot Island, 4 Disintegration, 5 Sea Ice, 6 Lataday Island, 7 Wilkins Ice Shelf, 8 Unidentified Regions, 9 Poor Fit");

%
%% PROBLEM 4) Optional: REDO CLASSIFICATION WITH DIFFERENT TRAINING SET
%

