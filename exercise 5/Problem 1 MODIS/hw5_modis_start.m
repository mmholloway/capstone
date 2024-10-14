%
% EPSc 407 Spring 2023
%
% This mfile provides skeleton code to perform  image classification
% of a 7-band MODIS scene.
%
clear, clf, clc;
%
% PROBLEM 1) READ IN THE IMAGE AND LOOK AT THE BANDS
% read in and look at all 7 bands with imagesc
%

red  = uint8(load('modis1.dat'));
nir  = uint8(load('modis2.dat')); % near infrared
blue  = uint8(load('modis3.dat'));
mir1 = uint8(load('modis4.dat')); % mid infrared
mir2 = uint8(load('modis5.dat'));
mir3  = uint8(load('modis6.dat'));
mir4 = uint8(load('modis7.dat'));

%
% look at each band separately and an rgb image
%

imshow(cat(3,red,mir1,blue));
imshow()

pause;
% CREATE A TRAINING SET
% Identify the location (row,column) of "training pixels" and assign them a "group number":
% You will need to identify at least 2 training pixels per group
%

% Each line of tpix should have 3 numbers in the following format:    
%     row,col,group
%      ^   ^  ^



%
% From these pixels, make a training set consisting of each training pixel's band values
% "train" should have the same number of rows as the number of training pixels, and the
% same number of columns as number of bands (in the MODIS case, 7).  
%   


% 
% PROBLEM 3A) PREPARE DATA FOR CLASSIFICATION
% Reshape image into one long vector of pixel band values.
% Convert from uint8 to double for classification.
%
   nx=400;
   ny=400;
   N=nx*ny;
   nz=7;
   
%
% PROBLEM 3B) CLASSIFY THE IMAGE
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

toc



%
% PROBLEM 3C) LOOK AT THE CLASSIFIED IMAGE
% Reshape the Class vector and Each group error vetor (Nx1) back into an nx x ny matrix.
%



%
% Visualize the classification.
% You could do this with any colormap, or you can make your own with RGB values
% for each Group Number that make sense to you. Feel free to use the map below
% or make your own.
%   
   map=[1,0.96863,0.92157;...          % Group 1: Clouds
        0.72941,0.83137,0.95686;...    % Group 2: Sea Ice
        0,0.74902,0.74902;...          % Group 3: Floating Ice Shelf
        0.95294,0.87059,0.73333;       % Group 4: Land Ice
        0.078431,0.16863,0.54902];     % Group 5: Open Ocean                    





%
% PROBLEM 3) OPTIONAL: LOOK AT THE MISFIT
%



%
% PROBLEM 4) Optional: REDO CLASSIFICATION WITH DIFFERENT TRAINING SET
%

