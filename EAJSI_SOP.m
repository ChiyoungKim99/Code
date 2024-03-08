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

load u200_201301_202312.mat
load lon.mat
load lat.mat

% remove global mean
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
lat_idx_3 = find(lat>5 & lat<=10);


jet_1 = squeeze(nansum(nansum(u200(lon_idx_1,lat_idx_1,:).*awm(lon_idx_1,lat_idx_1,:),1),2)./nansum(nansum(awm(lon_idx_1,lat_idx_1,:),1),2));
jet_2 = squeeze(nansum(nansum(u200(lon_idx_2,lat_idx_2,:).*awm(lon_idx_2,lat_idx_2,:),1),2)./nansum(nansum(awm(lon_idx_2,lat_idx_2,:),1),2));
jet_3 = squeeze(nansum(nansum(u200(lon_idx_3,lat_idx_3,:).*awm(lon_idx_3,lat_idx_3,:),1),2)./nansum(nansum(awm(lon_idx_3,lat_idx_3,:),1),2));

jet_idx = ((jet_1 - jet_2) + (jet_1 - jet_3))/2;


% djf_jet_idx = [];

% for i = 2013:2023
%     dec_idx = (i - 2013) * 12 + 12;
%     jan_idx = (i - 2013) * 12 + 1;
%     feb_idx = (i - 2013) * 12 + 2;
%     djf_jet_idx(i - 2012) = (jet_idx(dec_idx) + jet_idx(jan_idx) + jet_idx(feb_idx)) / 3;
% end

load lat_grid.mat;
load lon_grid.mat;

UTLS_path = '/data01/satellite/IASI/O3/PROFILE/METOPB/UTLS/';
Strato_path = '/data01/satellite/IASI/O3/PROFILE/METOPB/Strato/';
% TCO_path = '/data01/satellite/IASI/O3/TCO/METOPC/REGRID_0.5/';
% TROPOMI_path = '/data01/satellite/TROPOMI/O3/TCO/REGRID_0.5/';

for year = 2014:2023
    for month = 1:12
        start_date = datenum([num2str(year) num2str(month,'%02d') '01'], 'yyyymmdd');
        end_date = eomdate(start_date); % Get the last day of the month
        date_range = start_date:end_date;

        sop = NaN(size(lon_grid, 1), size(lat_grid, 2), numel(date_range));

        for t = 1:numel(date_range)
            date_str = datestr(date_range(t), 'yyyymmdd');

            try
                load([UTLS_path 'IASI_REGRID_0.5_UTLS_' date_str '.mat']);
                utlso_ia_re(utlso_ia_re == 0) = NaN;
                load([Strato_path 'IASI_REGRID_0.5_Strato_' date_str '.mat']);
                so_ia_re(so_ia_re == 0) = NaN;
                sop(:, :, t) = utlso_ia_re ./ so_ia_re;
            catch
                disp(['No data for IASI TCO ' date_str])
            end
        end

        lat_asia = find(lat_grid(:,1) >= 30 & lat_grid(:,1) <= 55);
        lon_asia = find(lon_grid(1,:) >= 110 & lon_grid(1,:) <= 160);

        sop_asia = sop(lat_asia, lon_asia, :);

        idx = (sop_asia >= 0.25);

        value_idx = isfinite(sop_asia);

        range = nansum(value_idx, 3);

        freq(:,:, (year-2014)*12 + month) = 100 * nansum(idx, 3) ./ range;
    end
end

% nlat=length(lat_asia); nlon=length(lon_asia);

% freq = freq_month - reshape(repmat(squeeze(nanmean(reshape(freq_month,[nlat nlon 12 10]),4)),[1 1 1 10]),[nlat nlon 120]);


% for i = 2014:2023
% 	for j = 1:5
% 		jet_JFMAM(5 * (i-2014) + j) = jet_idx(12*(i - 2013) + j);
% 	end
% end



% for i = 2014:2023
% 	for j = 1:5
% 		freq_JFMAM(:,:,5 * (i-2014) + j) = freq(:,:,12*(i - 2014) + j);
% 	end
% end

for i = 2014:2023
    for j = 3:5
        jet_MAM(3 * (i-2014) + j - 2) = jet_idx(12 * (i - 2013) + j);
        freq_MAM(:,:,3 * (i-2014) + j - 2) = freq(:,:,12 * (i - 2014) + j);
    end
end





nlat=length(lat_asia); nlon=length(lon_asia);
temp1 = zscore(jet_MAM); temp2 = freq_MAM;
temp1(isnan(temp1))=0; temp2(isnan(temp2))=0;
for i = 1:nlon
    for j = 1:nlat
        tval = regstats(squeeze(temp2(j,i,:)), temp1(:),'linear',{'tstat'} ); 
        regcf(j,i) = tval.tstat.beta(2); 
        regcf_pval(j,i) = tval.tstat.pval(2); 
    end
end



load coast_cykim.mat

figure('Position',[50 50 1400 600]); h=gca;
hold on;
% contourf(lon_grid(lat_asia, lon_asia), lat_grid(lat_asia, lon_asia), regcf, 'linestyle', 'none');
% p = pcolor(lon_grid(lat_asia, lon_asia), lat_grid(lat_asia, lon_asia), regcf);
% set(p, 'EdgeColor', 'none');
imagescn(lon_grid(lat_asia, lon_asia), lat_grid(lat_asia, lon_asia), regcf);
% colormap(french(20))
colormap(nclCM(322,20))
a = colorbar('YTick', -10:2:10);
% a.Label.String = 'SOP Frequency / EAJSI [%]';
a.Label.Rotation = 90;  % Set the rotation angle to 0 degrees for the colorbar label
a.FontSize = 16;
caxis([-10 10])
mask = regcf_pval <= 0.05;
stipple(lon_grid(lat_asia, lon_asia), lat_grid(lat_asia, lon_asia), mask, 'density',200)
plot(lonmap, latmap, 'k', 'LineWidth', 2);
xlim([110 160]); ylim([30 55]);
h.XTick = 60:15:170; h.XTickLabel = {'60\circE','75\circE','90\circE','105\circE','120\circE','135\circE','150\circE','165\circE'};
h.YTick = 30:5:55; h.YTickLabel = {'30\circN','35\circN','40\circN','45\circN','50\circN','55\circN'};
h.XTickLabelRotation = 0; 
h.FontSize=16; 
box on; grid on;
title('SOP Frequency MAM','FontSize',24, 'FontWeight', 'bold');

saveas(gcf, ['/home/cykim/figures/SOP_Frequency_regression_MAM.png']);




