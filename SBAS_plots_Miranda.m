%close all

%plot deformation time series for each site

ticks = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ...
    23 24 25 26 27 28 29];
%ticks = [1 4 6 9 11 13 16 18 21 24 26 29]
labels = {'06/16/2020', '06/28/2020', '07/10/2020', '07/22/2020', '08/03/2020'...
    '08/15/2020', '08/27/2020', '09/08/2020', '09/20/2020', '06/11/2021'...
    '06/23/2021', '07/05/2021', '07/17/2021', '07/29/2021', '08/10/2021'...
    '08/22/2021', '09/03/2021', '09/15/2021', '09/27/2021', '06/06/2022'...
    '06/18/2022', '06/30/2022', '07/12/2022', '07/24/2022', '08/05/2022'...
    '08/17/2022', '08/29/2022', '09/10/2022', '09/22/2022'};

% Below is the structure for plotting the deformation time series for a
% given pixel. Deformation is relative to the reference pixels chosen. We
% chose reference pixels that display high coherence over the entire
% 51-interferogram set. These reference pixels are located near
% stream/river and thermokarst lake banks and are interpreted to be
% gravel/sandy deposits that don't experience subsidence during the summer.

%Figures below are for pixels referenced to ice core sites, but can be
%changed to any pixel of interest for the
%"plot(squeeze(time_series(x,x:)))" lines
%% Set of 3-rowed time series

figure %T1-CHRONO-CPC-1
tiledlayout(3,1)
nexttile
plot(squeeze(time_series(346,403,:)))
title('Relative Elevation for T1-CHRONO-CPC-1')
ylabel('Relative Elevation (cm)')
xlabel('Time')
set(gca,"XTick",ticks)
set(gca,'XTickLabels', labels)
ylim([-10 5])

%figure %T2-CHRONO-CPC-2 ICE-POOR
nexttile
plot(squeeze(time_series(514,386,:)))
title('Relative Elevation for T2-CHRONO-CPC-2')
ylabel('Relative Elevation (cm)')
xlabel('Time')
set(gca,"XTick",ticks)
set(gca,'XTickLabels', labels)
ylim([-10 5])

%figure %T3-CHRONO-CPC-6 ICE-POOR
nexttile
plot(squeeze(time_series(418,360,:)))
title('Relative Elevation for T3-CHRONO-CPC-6')
ylabel('Relative Elevation (cm)')
xlabel('Time')
set(gca,"XTick",ticks)
set(gca,'XTickLabels', labels)
ylim([-10 5])

figure
tiledlayout(3,1)
%figure %T4-CHRONO-CPC-3 ICE-POOR
nexttile
plot(squeeze(time_series(192,298,:)))
title('Relative Elevation for T4-CHRONO-CPC-3')
ylabel('Relative Elevation (cm)')
xlabel('Time')
set(gca,"XTick",ticks)
set(gca,'XTickLabels', labels)
ylim([-10 5])

%figure %T5-CHRONO-CPC-7 ICE-RICH
nexttile
plot(squeeze(time_series(235,370,:)))
title('Relative Elevation for T5-CHRONO-CPC-7')
ylabel('Relative Elevation (cm)')
xlabel('Time')
set(gca,"XTick",ticks)
set(gca,'XTickLabels', labels)
ylim([-10 5])

%figure %T6-CHRONO-CPC-8 ICE-RICH
nexttile
plot(squeeze(time_series(264,374,:)))
title('Relative Elevation for T6-CHRONO-CPC-8')
ylabel('Relative Elevation (cm)')
xlabel('Time')
set(gca,"XTick",ticks)
set(gca,'XTickLabels', labels)
ylim([-10 5])

figure
tiledlayout(3,1)
%figure %T8-DREW Point ICE-RICH
nexttile
plot(squeeze(time_series(276,161,:)))
title('Relative Elevation for T8-Drew Point')
ylabel('Relative Elevation (cm)')
xlabel('Time')
set(gca,"XTick",ticks)
set(gca,'XTickLabels', labels)
ylim([-10 5])

%figure %T9-DREW Point ICE-RICH
nexttile
plot(squeeze(time_series(275,161,:)))
title('Relative Elevation for T9-Drew Point')
ylabel('Relative Elevation (cm)')
xlabel('Time')
set(gca,"XTick",ticks)
set(gca,'XTickLabels', labels)
ylim([-10 5])

%figure %T10-East Teshekpuk ICE-POOR
nexttile
plot(squeeze(time_series(1123,694,:)))
title('Relative Elevation for T10-East Teshepuk')
ylabel('Relative Elevation (cm)')
xlabel('Time')
set(gca,"XTick",ticks)
set(gca,'XTickLabels', labels)
ylim([-10 5])

%%%%%%%%%%%%%%%%%%%%
%% Individual 3x3 grid for viewing
figure %T1-CHRONO-CPC-1
tiledlayout(3,3)
nexttile
plot(squeeze(time_series(346,403,:)))
title('Relative Elevation for T1-CHRONO-CPC-1')
ylabel('Relative Elevation (cm)')
xlabel('Time')
set(gca,"XTick",ticks)
set(gca,'XTickLabels', labels)
ylim([-10 5])

%figure %T2-CHRONO-CPC-2 ICE-POOR
nexttile
plot(squeeze(time_series(514,386,:)))
title('Relative Elevation for T2-CHRONO-CPC-2')
ylabel('Relative Elevation (cm)')
xlabel('Time')
set(gca,"XTick",ticks)
set(gca,'XTickLabels', labels)
ylim([-10 5])

%figure %T3-CHRONO-CPC-6 ICE-POOR
nexttile
plot(squeeze(time_series(418,360,:)))
title('Relative Elevation for T3-CHRONO-CPC-6')
ylabel('Relative Elevation (cm)')
xlabel('Time')
set(gca,"XTick",ticks)
set(gca,'XTickLabels', labels)
ylim([-10 5])

%figure %T4-CHRONO-CPC-3 ICE-POOR
nexttile
plot(squeeze(time_series(192,298,:)))
title('Relative Elevation for T4-CHRONO-CPC-3')
ylabel('Relative Elevation (cm)')
xlabel('Time')
set(gca,"XTick",ticks)
set(gca,'XTickLabels', labels)
ylim([-10 5])

%figure %T5-CHRONO-CPC-7 ICE-RICH
nexttile
plot(squeeze(time_series(235,370,:)))
title('Relative Elevation for T5-CHRONO-CPC-7')
ylabel('Relative Elevation (cm)')
xlabel('Time')
set(gca,"XTick",ticks)
set(gca,'XTickLabels', labels)
ylim([-10 5])

%figure %T6-CHRONO-CPC-8 ICE-RICH
nexttile
plot(squeeze(time_series(264,374,:)))
title('Relative Elevation for T6-CHRONO-CPC-8')
ylabel('Relative Elevation (cm)')
xlabel('Time')
set(gca,"XTick",ticks)
set(gca,'XTickLabels', labels)
ylim([-10 5])

%figure %T8-DREW Point ICE-RICH
nexttile
plot(squeeze(time_series(276,161,:)))
title('Relative Elevation for T8-Drew Point')
ylabel('Relative Elevation (cm)')
xlabel('Time')
set(gca,"XTick",ticks)
set(gca,'XTickLabels', labels)
ylim([-10 5])

%figure %T9-DREW Point ICE-RICH
nexttile
plot(squeeze(time_series(275,161,:)))
title('Relative Elevation for T9-Drew Point')
ylabel('Relative Elevation (cm)')
xlabel('Time')
set(gca,"XTick",ticks)
set(gca,'XTickLabels', labels)
ylim([-10 5])

%figure %T10-East Teshekpuk ICE-POOR
nexttile
plot(squeeze(time_series(1123,694,:)))
title('Relative Elevation for T10-East Teshepuk')
ylabel('Relative Elevation (cm)')
xlabel('Time')
set(gca,"XTick",ticks)
set(gca,'XTickLabels', labels)
ylim([-10 5])


%% Combining all time series together

figure %combined
plot(squeeze(time_series(346,403,:)),'r')
hold on
plot(squeeze(time_series(514,386,:)),'g')
hold on
plot(squeeze(time_series(418,360,:)),'b')
hold on
plot(squeeze(time_series(192,298,:)),'c')
hold on
plot(squeeze(time_series(235,370,:)),'m')
hold on
plot(squeeze(time_series(264,374,:)),'y')
hold on
orange = [1 0.5 0];
plot(squeeze(time_series(276,161,:)),'Color',orange)
hold on
plot(squeeze(time_series(275,161,:)),'w')
hold on
gray = [.7 .7 .7];
plot(squeeze(time_series(1123,694,:)),'Color',gray)
title('Relative Elevation for Each Ice Core Site')
set(gca,"XTick",ticks)
set(gca,'XTickLabels', labels)
legend({'T1-CHRONO-CPC-1 ICE-POOR','T2-CHRONO-CPC-2 ICE-POOR','T3-CHRONO-CPC-6 ICE-POOR',...
    'T4-CHRONO-CPC-3 ICE-POOR','T5-CHRONO-CPC-7 ICE-RICH','T6-CHRONO-CPC-8 ICE-RICH',...
    'T8-DREW Point ICE-RICH','T9-DREW Point ICE-RICH','T10-East Teshekpuk ICE-POOR'},...
    'Location','northwest')
ylim([-10 5])

%% Reference pixels checking plot

figure
imAlpha=ones(size(mask));
imAlpha(isnan(mask))=0;
imagesc(avecc','AlphaData',imAlpha');
set(gca,'color',0*[1 1 1]);
title('Reference Pixels Check')
c = colorbar;
xlabel(c, 'Coherence [-]')
set(gca,'XTick',[378 756 1134 1512 1890])      %~pixel values corresponding to long.
set(gca,'XTickLabel',[-153.79 -153.37 -152.95 -152.53 -152.11])
set(gca,'YTick',[254 508 762 1016])      %~pixel values corresponding to lat.
set(gca,'YTickLabel',[70.81 70.67 70.53 70.38])
xlabel('Longitude')
ylabel('Latitude')
axis on
hold on
% ref pixel locations below
plot(1241, 597, 'rx','Linewidth', 2,'MarkerSize',10)
plot(1185, 186, 'rx','Linewidth', 2,'MarkerSize',10)
plot(183, 923, 'rx','Linewidth', 2,'MarkerSize',10)
plot(857, 859, 'rx','Linewidth', 2,'MarkerSize',10)

