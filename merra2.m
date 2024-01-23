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

MERRA_path = '/data01/MERRA2/monthly/';

start_date = datetime('20130101', 'Format', 'yyyyMMdd');
end_date = datetime('20231231', 'Format', 'yyyyMMdd');
date_range = start_date : calmonths(1) : end_date;


u200 = NaN(576, 361, numel(date_range));

for t = 1:numel(date_range)
    tic
    date_str = datestr(date_range(t), 'yyyymm');
    try
        u = ncread([MERRA_path 'MERRA2_400.instM_3d_asm_Np.' date_str '.nc4'], 'U');
        u200(:, :, t) = squeeze(u(:,:,23));
        disp(date_str);
    catch
        u = ncread([MERRA_path 'MERRA2_401.instM_3d_asm_Np.' date_str '.nc4'], 'U');
        u200(:, :, t) = squeeze(u(:,:,23));
        disp(date_str);
    end 
    toc
end

save_file_path = '/data01/MERRA2/monthly/u200_201301_202312.mat';
save(save_file_path, 'u200');


lon = ncread('/data01/MERRA2/monthly/MERRA2_400.instM_3d_asm_Np.201301.nc4', 'lon');
lat = ncread('/data01/MERRA2/monthly/MERRA2_400.instM_3d_asm_Np.201301.nc4', 'lat');

lat_file_path = '/data01/MERRA2/monthly/lat.mat';
save(lat_file_path, 'lat');
lon_file_path = '/data01/MERRA2/monthly/lon.mat';
save(lon_file_path, 'lon');