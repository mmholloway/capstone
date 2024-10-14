%% Miranda's Edited Version of Roger's post_process_data_ascending.m

%% Load all of the InSAR files you are considering
clear all; close all; clc;

% Added by Miranda
subdir = 'subdir_70/subdir_1718';
addpath(strcat('C:\Users\mmpho\sent_test\',subdir))

% Info from dem.rsc
nr=1890; naz=1016; % image size:  nr = # of x (range) pixels (width);
% naz = # of y (azimuth) pixels (file_length)
n=68; % number of slcs
% Last semester's dataset was 1468 pixels for both x and y, 45 SLCs

%Just read in coherent inteferograms with a temporal
%baseline of 48 days or less
Bperp0=load('Bperp.out');
Tm0=load('Tm.out');
deltime0=load('deltime.out');
timedeltas0=load('timedeltas.out');


Bperp1=zeros(size(Bperp0));
Tm1=zeros(size(Tm0));
deltime1=zeros(size(deltime0));


cells2=importdata('sbas_list');
%cells2=importdata('sbas_list36_ascending');
N2=length(cells2);

sbas_dates0=string(cells2.textdata);
sbas_data0=string(cells2.data);
cells=importdata('intlist');
%cells=importdata('intlist36');
N=length(cells);
lambda=5.6; %wavelength

unw_phase=zeros(nr,naz,N);
%uni=zeros(nr,naz,N);
amps=zeros(nr,naz,N);
ints=zeros(nr,naz,N);
coh=zeros(nr,naz,N);
%uni=zeros(nr,naz,N);
date_pair=cell(2,N);
doy_pair=cell(2,N);

%Read in the unwrapped phase (unw), coherence (coh), amplitude (amp) and
%unimodally-corrected unwrapped phase (uni)
for i=1:N
    disp(i)
    strint=cells{i};
    % strint1=strcat('ints/',strint);
    strint1=strint;
    strunw1=strrep(strint,'.int','.unw');
    % strunw=strcat('unws/',strunw1);
    strunw=strunw1;
    stramp1=strrep(strint,'.int','.amp');
    % stramp=strcat('amps/',stramp1);
    stramp=stramp1;
    strcc1=strrep(strint,'.int','.cc');
    % strcc=strcat('ccs/',strcc1);
    strcc=strcc1;
    struni1=strrep(strint,'.int','.int');
    % struni=strcat('ints/',struni1);
    struni=struni1;
    % correlations
    filename_c=sprintf('%s',strcc);
    fid=fopen(filename_c);
    dat=fread(fid,[2*nr,inf],'float','ieee-le');
    temp=dat((nr+1):end,:);
    coh(:,:,i)=temp;
    fclose(fid);
    % unwrapped phase
    filename=sprintf('%s',strunw);
    fid=fopen(filename);
    dat=fread(fid,[2*nr,inf],'float','ieee-le');
    temp=dat(nr+1:end,:);
    unw_phase(:,:,i)=temp;
    fclose(fid);
    % Interferograms
    %filename=sprintf('%s',struni);
    %fid=fopen(filename);
    %dat=fread(fid,[2*nr,inf],'float','ieee-le');
    %temp=dat(1:2:end,1:naz)+1i*dat(2:2:end,1:naz);
    %phase(:,:,i)=temp;
    %fclose(fid);
    % wrapped phase
    filename=sprintf('%s',strint1);
    fid=fopen(filename);
    dat=fread(fid,[2*nr,inf],'float','ieee-le');
    temp=dat(1:2:2*nr-1,:)+1i*dat(2:2:2*nr,:);
    ints(:,:,i)=temp;
    fclose(fid);
    % Amplitude
    filename=sprintf('%s',stramp);
    fid=fopen(filename);
    dat=fread(fid,[2*nr,inf],'float','ieee-le');
    temp=dat(1:2:2*nr-1,:)+1i*dat(2:2:2*nr,:);
    amps(:,:,i)=temp;
    fclose(fid);
    % date information
    split1=strsplit(strint,'_');
    strint2=split1{2};
    split2=strsplit(strint2,'.');
    d1=split1{1};
    d2=split2{1};

    date1=strcat(d1(5:6),'/',d1(7:8),'/',d1(1:4));
    date2=strcat(d2(5:6),'/',d2(7:8),'/',d2(1:4));

    date1_vec=datetime(date1,'InputFormat','MM/dd/yyyy');
    date2_vec=datetime(date2,'InputFormat','MM/dd/yyyy');

    doy1=day(date1_vec,'dayofyear');
    doy2=day(date2_vec,'dayofyear');

    date_pair{1,i}=date1;
    date_pair{2,i}=date2;
    doy_pair{1,i}=doy1;
    doy_pair{2,i}=doy2;
end

%% Get an average coherence file
avecc = mean(coh,3);
cor_mask=avecc;
alpha = 0.44; %0.17
cor_mask(cor_mask<alpha)=nan;
cor_mask(cor_mask>alpha)=1;

%generated coherence mask
mask=cor_mask;

%% Now we will do the atmospheric correction
%remove topographically-correlated atmospheric noise from interferograms
% Depending on how large the scene is, you may want to base the correction
% solely on a few areas with the topographic relief changes quite a bit.
% But if the scene is small, you can indeed use the whole scene.
% This bit also requires a calibration pixel (or set of pixels) to ensure
% that all of the scenes are set to the same "datum"

% Load the dem
nr0=14401; %%%%%%%%%%%%%%%%%%%%%%%% where do these values come from?? (something with DEM size?)
naz0=7201;
fid=fopen('elevation.dem','r');
dem0=fread(fid,[nr0,naz0],'int16'); % x length first, y length second
fclose(fid);

dem=imresize(dem0,[nr,naz]);
masked_unw_phase = unw_phase;%.*mask;
%masked_unw_phase1 = uni.*mask;

% Pick some pixels for calibration
pixels = [1140 1070]; % [range azimuth] %%%%%%%%%%%%%%%% choosing good values here? why?
[sz2,~] = size(pixels);

phase = zeros(nr,naz,N);
phase1 = zeros(nr,naz,N);
corrections=zeros(nr,naz,N);
phase_rshp = [];
dem_rshp = [];
for int = 1:N
    block_phase = masked_unw_phase(:,:,int);
    indx = isnan(block_phase);
    line = polyfit(dem(~indx),block_phase(~indx),1);
    correction = (line(1)*dem + line(2));
    corrections(:,:,int)=correction;
    %phase(:,:,int) = masked_unw_phase(:,:,int);
    phase(:,:,int) = masked_unw_phase(:,:,int) - correction;
end

% clear masked_unw_phase
% clear unw_phase
% clear ints

%% IG Error Plots code (from Alex thanks king)

%setting alpha channel to recognize NaN values in mask, setting to black
imAlpha=ones(size(mask));
imAlpha(isnan(mask))=0;
% set(gca,'color',0*[1 1 1]);

%looping through all interferograms
% CAREFUL, time taken will depend on # of ints and processing power
addpath('C:\Users\mmpho\OneDrive\Documents\radar-lab\sp24')
no_x_tix = 5;
no_y_tix = 5;
[pix, labels] = igram_latlong(subdir,no_x_tix,no_y_tix);
no_x_tix = no_x_tix + 1;
no_y_tix = no_y_tix + 1;

if (~isfolder('C:\Users\mmpho\OneDrive - Washington University in St. Louis\Year 4\Capstone\IG Error Plots'))
    mkdir('C:\Users\mmpho\OneDrive - Washington University in St. Louis\Year 4\Capstone\IG Error Plots');
end

for i= 1:N % [17 35 44]
    fig = figure('Name',strcat("Alpha value = ",num2str(alpha)));
    tiledlayout(2,2)
    nexttile
    % INT FIG w/ MASK %
    imagesc(angle(ints(:,:,i).*mask).','AlphaData',imAlpha')
    title(['Wrapped Phase #',num2str(i)],[date_pair{1,i},' - ', date_pair{2,i}])
    colorbar
    % probably a better way to do this, but setting manually lat/lon
    % will need to change this depending on your region
    set(gca,'XTick',pix(1,1:no_x_tix))      %~pixel values corresponding to long.
    set(gca,'XTickLabel',labels(1,1:no_x_tix))
    set(gca,'YTick',pix(2,1:no_y_tix))      %~pixel values corresponding to lat.
    set(gca,'YTickLabel',labels(2,1:no_y_tix))
    xlabel('Longitude (degrees)')
    ylabel('Latitude (degrees)')
    set(gca,'FontSize',16)

    % COH FIG %
    nexttile
    imagesc(coh(:,:,i).')
    title(['Coherence #',num2str(i)],[date_pair{1,i},' - ', date_pair{2,i}])
    colorbar
    % probably a better way to do this, but setting manually lat/lon
    % will need to change this depending on your region
    set(gca,'XTick',pix(1,1:no_x_tix))      %~pixel values corresponding to long.
    set(gca,'XTickLabel',labels(1,1:no_x_tix))
    set(gca,'YTick',pix(2,1:no_y_tix))      %~pixel values corresponding to lat.
    set(gca,'YTickLabel',labels(2,1:no_y_tix))
    xlabel('Longitude (degrees)')
    ylabel('Latitude (degrees)')
    set(gca,'FontSize',16)

    % UNW PHASE w/ MASK %
    nexttile
    imagesc((unw_phase(:,:,i).*mask.*lambda/(4*pi)).','AlphaData',imAlpha')
    title(['Unwrapped Phase #',num2str(i)],[date_pair{1,i},' - ', date_pair{2,i}])
    colorbar
    xlabel(colorbar, '[cm]')
    % probably a better way to do this, but setting manually lat/lon
    % will need to change this depending on your region
    set(gca,'XTick',pix(1,1:no_x_tix))      %~pixel values corresponding to long.
    set(gca,'XTickLabel',labels(1,1:no_x_tix))
    set(gca,'YTick',pix(2,1:no_y_tix))      %~pixel values corresponding to lat.
    set(gca,'YTickLabel',labels(2,1:no_y_tix))
    xlabel('Longitude (degrees)')
    ylabel('Latitude (degrees)')
    set(gca,'FontSize',16)

    % % this is setting limits on the color bar (helpful to bring bounds in)
    stats_array = squeeze(unw_phase(:,:,i).*mask);
    row_stats_array = (stats_array(:));

    percent5 = prctile(row_stats_array,5);
    percent95 = prctile(row_stats_array,95);
    clim([percent5 percent95])

    % UNW PHASE HISTOGRAM DISTRIBUTION %
    nexttile
    histogram((unw_phase(:,:,i).*mask).')
    title(['Unwrapped Phase Distribution #',num2str(i)],[date_pair{1,i},' - ', date_pair{2,i}])
    set(gcf, 'Position', get(0, 'Screensize')); %forcing  automatic full screen resolution
    %figure

    % Save file
    exportgraphics(fig,strcat('C:\Users\mmpho\OneDrive - Washington University in St. Louis\Year 4\Capstone\IG Error Plots\',num2str(alpha),'_',num2str(i),'.png'))
    close;
end

%% Plot some unwrapped phase distributions for raw images

% Good images: 1, 3, 4, 5, 6, 9
% Bad images: 21, 22, 28, 29, 31, 32, 35-40 (very low coherence), 41, 43,
% 45, 46, 49, 50, 51

for i = [1, 3, 4, 5, 6, 9]
    % Plot good image UW phase distributions
    fig = figure;
    tiledlayout(1,2)
    % nexttile
    % histogram((unw_phase(:,:,i)).')
    % title(['Good, No Mask UW Phase Distribution #',num2str(i)],[date_pair{1,i},' - ', date_pair{2,i}])
    % xlabel("Unwrapped Phase (radians)")
    % ylabel("Number of Samples")
    % set(gcf, 'Position', get(0, 'Screensize'));
    % axis tight
    % 
    % nexttile
    % histogram((unw_phase(:,:,i).*mask).')
    % title(['Good, Masked UW Phase Distribution #',num2str(i)],[date_pair{1,i},' - ', date_pair{2,i}])
    % xlabel("Unwrapped Phase (radians)")
    % ylabel("Number of Samples")
    % set(gcf, 'Position', get(0, 'Screensize'));
    % axis tight

    nexttile
    histogram((unw_phase(:,:,i).*lambda/(4*pi)).')
    title(['Good, No Mask UW Phase Distribution #',num2str(i)],[date_pair{1,i},' - ', date_pair{2,i}])
    xlabel("Unwrapped Phase (cm)")
    ylabel("Number of Samples")
    set(gcf, 'Position', get(0, 'Screensize'));
    axis tight

    nexttile
    histogram((unw_phase(:,:,i).*mask.*lambda/(4*pi)).')
    title(['Good, Masked UW Phase Distribution #',num2str(i)],[date_pair{1,i},' - ', date_pair{2,i}])
    xlabel("Unwrapped Phase (cm)")
    ylabel("Number of Samples")
    set(gcf, 'Position', get(0, 'Screensize'));
    axis tight

    exportgraphics(fig,strcat('C:\Users\mmpho\OneDrive - Washington University in St. Louis\Year 4\Capstone\IG Error Plots\','Good_WWO_Mask_',num2str(i),'.png'))
    close;
end

for i = [21, 22, 28, 29, 31, 32]
    % Plot bad image UW phase distributions
    fig = figure;
    tiledlayout(1,2)
    % nexttile
    % histogram((unw_phase(:,:,i)).')
    % title(['Bad, No Mask UW Phase Distribution #',num2str(i)],[date_pair{1,i},' - ', date_pair{2,i}])
    % xlabel("Unwrapped Phase (radians)")
    % ylabel("Number of Samples")
    % set(gcf, 'Position', get(0, 'Screensize'));
    % axis tight
    % 
    % nexttile
    % histogram((unw_phase(:,:,i).*mask).')
    % title(['Bad, Masked UW Phase Distribution #',num2str(i)],[date_pair{1,i},' - ', date_pair{2,i}])
    % xlabel("Unwrapped Phase (radians)")
    % ylabel("Number of Samples")
    % set(gcf, 'Position', get(0, 'Screensize'));
    % axis tight

    nexttile
    histogram((unw_phase(:,:,i).*lambda/(4*pi)).')
    title(['Bad, No Mask UW Phase Distribution #',num2str(i)],[date_pair{1,i},' - ', date_pair{2,i}])
    xlabel("Unwrapped Phase (cm)")
    ylabel("Number of Samples")
    set(gcf, 'Position', get(0, 'Screensize'));
    axis tight

    nexttile
    histogram((unw_phase(:,:,i).*mask.*lambda/(4*pi)).')
    title(['Bad, Masked UW Phase Distribution #',num2str(i)],[date_pair{1,i},' - ', date_pair{2,i}])
    xlabel("Unwrapped Phase (cm)")
    ylabel("Number of Samples")
    set(gcf, 'Position', get(0, 'Screensize'));
    axis tight

    exportgraphics(fig,strcat('C:\Users\mmpho\OneDrive - Washington University in St. Louis\Year 4\Capstone\IG Error Plots\','Bad_WWO_Mask_',num2str(i),'.png'))
    close;
end