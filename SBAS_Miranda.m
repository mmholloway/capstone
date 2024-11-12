%%Convention is Master_slave.unw, and crossmul takes master*conj(slave),
%%such that a positive values in the unwrapped phase corresponds to a
%%subsidence 


%% Load all of the InSAR files you are considering
clear all; clc; % close all;

% Added by Miranda
subdir = 'sbas_24';
addpath(strcat('C:\Users\mmpho\sent_test\',subdir))

nr=1890; naz=1016;% image size:  nr = # of x (range) pixels; 
                             % naz = # of y (azimuth) pixels
cells=importdata('intlist');
N=length(cells);
lambda=5.6; %wavelength

unw_phase=zeros(nr,naz,N);
ints=zeros(nr,naz,N);
coh=zeros(nr,naz,N);
date_pair=cell(2,N);
doy_pair=cell(2,N);

%Read in the unwrapped phase (unw), coherence (coh), and unwrapped phase (int)
for i=1:N
    disp(i);
    strint1=cells{i};
    %strint=strcat(strint1,'/filt_fine.int');
    strcor=strrep(strint1,'.int','.cc');
    strunw=strrep(strint1,'.int','.unw');
    
% % correlations
    filename_c=sprintf('%s',strcor); 
    fid=fopen(filename_c);
    dat=fread(fid,[2*nr,inf],'float','ieee-le');
    temp=dat((nr+1):end,:);
    %temp=dat;
    coh(:,:,i)=temp;
    fclose(fid);
%unwrapped phase
    filename=sprintf('%s',strunw);  
    fid=fopen(filename);
    dat=fread(fid,[2*nr,naz],'float','ieee-le');
    %temp1=dat(:,1:naz);
    temp2=dat((nr+1):end,:);
    unw_phase(:,:,i)=temp2;
    fclose(fid);
%wrapped phase
    filename=sprintf('%s',strint1);  
    fid=fopen(filename);
    dat=fread(fid,[2*nr,naz],'float','ieee-le');
    %temp1=dat(:,1:naz);
    temp2=dat((nr+1):end,:);
    ints(:,:,i)=dat(1:2:end,:)+1i*dat(2:2:end,:);
    fclose(fid);
% date information
split1=strsplit(strint1,'_');
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

% Get an average coherence file 
avecc = mean(coh,3);
cor_mask=avecc;
cor_mask(cor_mask<0.5)=nan;
cor_mask(cor_mask>0.5)=1;
mask=cor_mask;

%% Read in SBAS parameters and make them consistent with the other data
Bperp0=load('Bperp.out');
Tm0=load('Tm.out');
deltime0=load('deltime.out'); 
timedeltas0=load('timedeltas.out');

cells2=importdata('sbas_list');
N2=length(cells2);

sbas_dates0=string(cells2.textdata);
sbas_data0=string(cells2.data);

Bperp1=Bperp0;
Tm1=Tm0;
deltime1=deltime0;
timedeltas1=timedeltas0;

cellsn=importdata('geolist');
n=length(cellsn);

%get dayofyear information for each SLC
for i=1:n
    strint=cellsn{i};
    d1=strint(21:28);
    date1=strcat(d1(5:6),'/',d1(7:8),'/',d1(1:4));
    date1_v=datetime(date1,'InputFormat','MM/dd/yyyy');
    doy1_v=day(date1_v,'dayofyear');
    geo_date{i}=date1_v;
    geo_doy{i}=doy1_v;
end

%% Calculate average rate
%Use four pixels with high average coherence as phase reference pixels;
%this can be changed in the future if needed
r_ref=[1241 1185 183 857];   %  reference pixel location
az_ref=[597 186 923 859];  %  reference pixel location


phase=unw_phase;
rates=zeros(nr,naz,N);
refphase=mean(mean(phase(r_ref,az_ref,:),1),2);
%write out average rate for inteferograms of temporal baseline X or less
for i=1:N
    rates(:,:,i)=(phase(:,:,i)-refphase(i)).*(-lambda/(4*pi))./deltime1(i,2);
end
rate=mean(rates,3);

%convert mean deformation rate to units of cm/year
rate_cmyr=rate.*365;


%% Do True SBAS, then stack the velocities

Bperp=Bperp1;
Tm=Tm1;
deltime=deltime1;
velocities=zeros(nr,naz,size(Tm,2));

Tmi=pinv(Tm);

for i=1:nr
    disp(i);
    for j=1:naz
        pluck=squeeze(phase(i,j,:))-squeeze(refphase);
        sol=Tmi*pluck;
        velocities(i,j,:)=sol;
        
    end
end
        

phi_vec=zeros(nr,naz,length(geo_doy));
for i=1:nr
    disp(i);
    for j=1:naz
        phi=zeros(length(timedeltas0),1);
        for kk=2:n
            phi(kk)=phi(kk-1)+timedeltas0(kk-1).*velocities(i,j,kk-1);
        end
        phi_vec(i,j,:)=phi;
    end
end

time_series=phi_vec.*(-lambda/(4*pi)); %this is the vector displaying
% time series of deformation.

% Note that deformation reported by SBAS is relative to the chosen reference
% pixels, which we assume to have steady elevation throughout the season
% due to high coherence and a prior information of little subsidence
% occurring at those pixels.