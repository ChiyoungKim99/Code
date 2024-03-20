close all; clear all; clc

addpath /home/cykim/matlab/cdt
addpath /home/cykim/matlab/cdt/cdt_data
addpath /home/cykim/matlab/cdt/doc
addpath /home/cykim/matlab/func/mexcdf/mexnc
addpath /home/cykim/matlab/func/mexcdf/snctools
addpath /home/cykim/matlab/func
addpath /home/cykim/matlab/func/landmask
addpath /home/cykim/matlab/func/colormap
addpath /home/cykim/matlab/func/hatch
addpath /home/cykim/matlab/func/nclCM
addpath /home/cykim/grid
addpath /data01/MERRA2/monthly

MERRA_path = '/data01/MERRA2/monthly/';

start_date = datetime('20140101', 'Format', 'yyyyMMdd');
end_date = datetime('20231231', 'Format', 'yyyyMMdd');
date_range = start_date : calmonths(1) : end_date;

u200 = NaN(576, 361, numel(date_range));
v200 = NaN(576, 361, numel(date_range));



for t = 1:numel(date_range)
    tic
    date_str = datestr(date_range(t), 'yyyymm');
    try
        u = ncread([MERRA_path 'MERRA2_400.instM_3d_asm_Np.' date_str '.nc4'], 'U');
        v = ncread([MERRA_path 'MERRA2_400.instM_3d_asm_Np.' date_str '.nc4'], 'V');
        u200(:, :, t) = squeeze(u(:,:,23));
        v200(:, :, t) = squeeze(v(:,:,23));
        disp(date_str);
    catch
        u = ncread([MERRA_path 'MERRA2_401.instM_3d_asm_Np.' date_str '.nc4'], 'U');
        v = ncread([MERRA_path 'MERRA2_401.instM_3d_asm_Np.' date_str '.nc4'], 'V');
        u200(:, :, t) = squeeze(u(:,:,23));
        v200(:, :, t) = squeeze(v(:,:,23));
        disp(date_str);
    end 
    toc
end


load u200_201301_202312.mat
load lon.mat
load lat.mat
load lev.mat


a = 6371e3; nlat = length(lat); nlon = length(lon); 
dx = a * repmat(cosd(lat),[1 nlon]) .* deg2rad(repmat(gradient(lon'),[nlat 1]));
dy = a * deg2rad(repmat(gradient(lat),[1 nlon]));
aw = dx.*dy;
awm = permute(repmat(aw,[1 1 size(u200,3)]),[2 1 3]);

lon_idx_1 = find(lon>100 & lon<=180);
lat_idx_1 = find(lat>30 & lat<=35);

lon_idx_2 = find(lon>70 & lon<=170);
lat_idx_2 = find(lat>50 & lat<=60);

lon_idx_3 = find(lon>90 & lon<=160);
lat_idx_3 = find(lat>-5 & lat<=10);


jet_1 = squeeze(nansum(nansum(u200(lon_idx_1,lat_idx_1,:).*awm(lon_idx_1,lat_idx_1,:),1),2)./nansum(nansum(awm(lon_idx_1,lat_idx_1,:),1),2));
jet_2 = squeeze(nansum(nansum(u200(lon_idx_2,lat_idx_2,:).*awm(lon_idx_2,lat_idx_2,:),1),2)./nansum(nansum(awm(lon_idx_2,lat_idx_2,:),1),2));
jet_3 = squeeze(nansum(nansum(u200(lon_idx_3,lat_idx_3,:).*awm(lon_idx_3,lat_idx_3,:),1),2)./nansum(nansum(awm(lon_idx_3,lat_idx_3,:),1),2));

jet_idx = ((jet_1 - jet_2) + (jet_1 - jet_3))/2;

load lat_grid.mat;
load lon_grid.mat;

lat_asia = find(lat >= 20 & lat <= 55);
lon_asia = find(lon >= 110 & lon <= 160);

nlat=length(lat_asia); nlon=length(lon_asia);

u200_MAM = NaN(576, 361, 30);
v00_MAM = NaN(576, 361, 30);
jet_MAM = NaN(30,1);

for i = 2014:2023
    for j = 3:5
        jet_MAM(3 * (i-2014) + j - 2) = jet_idx(12 * (i - 2013) + j);
        u200_MAM(:,:,3 * (i-2014) + j - 2) = u200(:,:,12 * (i - 2014) + j);
        v200_MAM(:,:,3 * (i-2014) + j - 2) = v200(:,:,12 * (i - 2014) + j);
    end
end

jet_index = zscore(jet_MAM);

wind_positive = sqrt(squeeze(nanmean(u200_MAM(:,:,(jet_index>0)),3))'.^2 + squeeze(nanmean(v200_MAM(:,:,(jet_index>0)),3))'.^2);

%-------------------------------------------------------------------------------------------------------------------------------------------------
load coast_cykim.mat

figure('Position',[50 50 1600 800]); h=gca;
hold on;
contourf([lon; 180],lat, [wind_positive wind_positive(:,1)], 20, 'linestyle', 'none')
contourf(lon+360,lat, wind_positive, 20, 'linestyle', 'none')

% colormap(nclCM(438))
colormap(nclCM(38,20))

a = colorbar('YTick', 0:5:50);
a.Label.String = '200hPa vector wind (m/s)';
a.Label.Rotation = 90;  % Set the rotation angle to 0 degrees for the colorbar label
a.FontSize = 14;
caxis([0 50])
xlim([80 220]); ylim([15 80]);
quiversc(lon,lat, nanmean(u200_MAM(:,:,(jet_index>0)),3)', nanmean(v200_MAM(:,:,(jet_index>0)),3)','k','density',100)
quiversc(lon+360,lat, nanmean(u200_MAM(:,:,(jet_index>0)),3)', nanmean(v200_MAM(:,:,(jet_index>0)),3)','k','density',100)
plot(lonmap, latmap, 'k', 'LineWidth', 2);
plot(lonmap+360, latmap, 'k', 'LineWidth', 2);
% xlim([0 360]); ylim([-90 90]);

h.XTick = 90:30:240; h.XTickLabel = {'90\circE','120\circE','150\circE','180\circE','150\circW'};
h.YTick = 0:10:80; h.YTickLabel = {'0\circ','10\circN','20\circN','30\circN','40\circN','50\circN','60\circN','70\circN','80\circN'};
h.XTickLabelRotation = 0; 
h.FontSize=16; 
box on; grid on;

saveas(gcf, ['/home/cykim/figures/Positive_wind_MAM_crop_200hPa.png']);
%-------------------------------------------------------------------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------------------------------------------------------------------
load coast_hjkim

figure('Position',[50 50 1600 800]); h=gca;
hold on;
contourf([lon; 180],lat, [wind_positive wind_positive(:,1)], 20, 'linestyle', 'none')
contourf(lon+360,lat, wind_positive, 20, 'linestyle', 'none')

% colormap(nclCM(438))
colormap(nclCM(38,20))

a = colorbar('YTick', 0:5:50);
a.Label.String = '200hPa vector wind (m/s)';
a.Label.Rotation = 90;  % Set the rotation angle to 0 degrees for the colorbar label
a.FontSize = 14;
caxis([0 50])
quiversc(lon,lat, nanmean(u200_MAM(:,:,(jet_index>0)),3)', nanmean(v200_MAM(:,:,(jet_index>0)),3)','k')
quiversc(lon+360,lat, nanmean(u200_MAM(:,:,(jet_index>0)),3)', nanmean(v200_MAM(:,:,(jet_index>0)),3)','k')
plot(lonmap, latmap, 'k', 'LineWidth', 2);
plot(lonmap+360, latmap, 'k', 'LineWidth', 2);
xlim([0 360]); ylim([-90 90]);
h.XTick = 0:60:360; h.XTickLabel = {'0\circ','60\circE','120\circE','180\circE','120\circW','60\circW','0\circ'};
h.YTick = -90:30:90; h.YTickLabel = {'90\circS','60\circS','30\circS','0\circ','30\circN','60\circN','90\circN'};
h.XTickLabelRotation = 0; 
h.FontSize=16; 
box on; grid on;

saveas(gcf, ['/home/cykim/figures/Positive_wind_MAM_200hPa.png']);


%-------------------------------------------------------------------------------------------------------------------------------------------------

wind_negative = sqrt(squeeze(nanmean(u200_MAM(:,:,(jet_index<0)),3))'.^2 + squeeze(nanmean(v200_MAM(:,:,(jet_index<0)),3))'.^2);

load coast_hjkim.mat

figure('Position',[50 50 1600 800]); h=gca;
hold on;
contourf([lon; 180],lat, [wind_negative wind_negative(:,1)], 20,'linestyle', 'none')
contourf(lon+360,lat, wind_negative, 20, 'linestyle', 'none')

% colormap(nclCM(438))
colormap(nclCM(38,20))
a = colorbar('YTick', 0:5:50);
a.Label.String = '200hPa vector wind (m/s)';
a.Label.Rotation = 90;  % Set the rotation angle to 0 degrees for the colorbar label
a.FontSize = 14;
caxis([0 50])

quiversc(lon,lat, nanmean(u200_MAM(:,:,(jet_index<0)),3)', nanmean(v200_MAM(:,:,(jet_index<0)),3)','k')
quiversc(lon+360,lat, nanmean(u200_MAM(:,:,(jet_index<0)),3)', nanmean(v200_MAM(:,:,(jet_index<0)),3)','k')
plot(lonmap, latmap, 'k', 'LineWidth', 2);
plot(lonmap+360, latmap, 'k', 'LineWidth', 2);
xlim([0 360]); ylim([-90 90]);
h.XTick = 0:60:360; h.XTickLabel = {'0\circ','60\circE','120\circE','180\circE','120\circW','60\circW','0\circ'};
h.YTick = -90:30:90; h.YTickLabel = {'90\circS','60\circS','30\circS','0\circ','30\circN','60\circN','90\circN'};
h.XTickLabelRotation = 0; 
h.FontSize=16; 
box on; grid on;

saveas(gcf, ['/home/cykim/figures/Negative_wind_MAM_200hPa.png']);

%-------------------------------------------------------------------------------------------------------------------------------------------------

load coast_cykim.mat

figure('Position',[50 50 1600 800]); h=gca;
hold on;
contourf([lon; 180],lat, [wind_negative wind_negative(:,1)], 20, 'linestyle', 'none')
contourf(lon+360,lat, wind_negative, 20, 'linestyle', 'none')

% colormap(nclCM(438))
colormap(nclCM(38,20))

a = colorbar('YTick', 0:5:50);
a.Label.String = '200hPa vector wind (m/s)';
a.Label.Rotation = 90;  % Set the rotation angle to 0 degrees for the colorbar label
a.FontSize = 14;
caxis([0 50])
xlim([80 220]); ylim([15 80]);
quiversc(lon,lat, nanmean(u200_MAM(:,:,(jet_index>0)),3)', nanmean(v200_MAM(:,:,(jet_index>0)),3)','k','density',100)
quiversc(lon+360,lat, nanmean(u200_MAM(:,:,(jet_index>0)),3)', nanmean(v200_MAM(:,:,(jet_index>0)),3)','k','density',100)
plot(lonmap, latmap, 'k', 'LineWidth', 2);
plot(lonmap+360, latmap, 'k', 'LineWidth', 2);
% xlim([0 360]); ylim([-90 90]);

h.XTick = 90:30:240; h.XTickLabel = {'90\circE','120\circE','150\circE','180\circE','150\circW'};
h.YTick = 0:10:80; h.YTickLabel = {'0\circ','10\circN','20\circN','30\circN','40\circN','50\circN','60\circN','70\circN','80\circN'};
h.XTickLabelRotation = 0; 
h.FontSize=16; 
box on; grid on;

saveas(gcf, ['/home/cykim/figures/Negative_wind_MAM_crop_200hPa.png']);


%-------------------------------------------------------------------------------------------------------------------------------------------------


wind_diff = wind_positive - wind_negative;

load coast_hjkim.mat

figure('Position',[50 50 1600 800]); h=gca;
hold on;
contourf([lon; 180],lat, [wind_diff wind_diff(:,1)], 32,'linestyle', 'none')
contourf(lon+360,lat, wind_diff, 32, 'linestyle', 'none')

% cmocean('balance')
colormap(nclCM(233,32))
a = colorbar('YTick', -15:1:15);
a.Label.String = '200hPa vector wind (m/s)';
a.Label.Rotation = 90;  % Set the rotation angle to 0 degrees for the colorbar label
a.FontSize = 14;
caxis([-16 16])

quiversc(lon,lat, nanmean(u200_MAM(:,:,(jet_index>0)),3)' - nanmean(u200_MAM(:,:,(jet_index<0)),3)', nanmean(v200_MAM(:,:,(jet_index>0)),3)' - nanmean(v200_MAM(:,:,(jet_index<0)),3)','k')
quiversc(lon+360,lat, nanmean(u200_MAM(:,:,(jet_index>0)),3)' - nanmean(u200_MAM(:,:,(jet_index<0)),3)', nanmean(v200_MAM(:,:,(jet_index>0)),3)' - nanmean(v200_MAM(:,:,(jet_index<0)),3)','k')
plot(lonmap, latmap, 'k', 'LineWidth', 2);
plot(lonmap+360, latmap, 'k', 'LineWidth', 2);
xlim([0 360]); ylim([-90 90]);
h.XTick = 0:60:360; h.XTickLabel = {'0\circ','60\circE','120\circE','180\circE','120\circW','60\circW','0\circ'};
h.YTick = -90:30:90; h.YTickLabel = {'90\circS','60\circS','30\circS','0\circ','30\circN','60\circN','90\circN'};
h.XTickLabelRotation = 0; 
h.FontSize=16; 
box on; grid on;

saveas(gcf, ['/home/cykim/figures/Difference_wind_MAM_200hPa.png']);


%-------------------------------------------------------------------------------------------------------------------------------------------------

load coast_cykim.mat

figure('Position',[50 50 1600 800]); h=gca;
hold on;
contourf([lon; 180],lat, [wind_diff wind_diff(:,1)], 32,'linestyle', 'none')
contourf(lon+360,lat, wind_diff, 32, 'linestyle', 'none')

% cmocean('balance')
colormap(nclCM(233,32))
a = colorbar('YTick', -15:1:15);
a.Label.String = '200hPa vector wind (m/s)';
a.Label.Rotation = 90;  % Set the rotation angle to 0 degrees for the colorbar label
a.FontSize = 14;
caxis([-16 16])
xlim([80 220]); ylim([15 80]);

quiversc(lon,lat, nanmean(u200_MAM(:,:,(jet_index>0)),3)' - nanmean(u200_MAM(:,:,(jet_index<0)),3)', nanmean(v200_MAM(:,:,(jet_index>0)),3)' - nanmean(v200_MAM(:,:,(jet_index<0)),3)','k','density',100)
quiversc(lon+360,lat, nanmean(u200_MAM(:,:,(jet_index>0)),3)' - nanmean(u200_MAM(:,:,(jet_index<0)),3)', nanmean(v200_MAM(:,:,(jet_index>0)),3)' - nanmean(v200_MAM(:,:,(jet_index<0)),3)','k','density',100)
plot(lonmap, latmap, 'k', 'LineWidth', 2);
plot(lonmap+360, latmap, 'k', 'LineWidth', 2);
% xlim([0 360]); ylim([-90 90]);

h.XTick = 90:30:240; h.XTickLabel = {'90\circE','120\circE','150\circE','180\circE','150\circW'};
h.YTick = 0:10:80; h.YTickLabel = {'0\circ','10\circN','20\circN','30\circN','40\circN','50\circN','60\circN','70\circN','80\circN'};
h.XTickLabelRotation = 0; 
h.FontSize=16; 
box on; grid on;

saveas(gcf, ['/home/cykim/figures/Difference_wind_MAM_crop_200hPa.png']);


