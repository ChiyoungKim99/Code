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
addpath /data01/satellite/IASI/O3/PROFILE/METOPB/ozonesonde_poh
addpath /data01/ozonesonde/pohang

sonde_path = '/data01/ozonesonde/pohang/';

r = 287.3;
g = 9.80665;
t0 = 273.15;
p0 = 1.01325e5;

altitude = linspace(1,41,41);
altitude_int = 0:0.1:41;

start_date = datenum('20130308', 'yyyymmdd');
end_date = datenum('20221231', 'yyyymmdd');
date_range = start_date:end_date;


for t = 1:numel(date_range)
    tic
    date_str = datestr(date_range(t), 'yyyymmdd');
	file_pattern = [sonde_path '*' date_str '*'];
	files = dir(file_pattern);

    try
        csv_data = readtable(files(1).name);
        load(['IASI_pohang_nearest_' date_str '.mat'])
        gph = csv_data.GPHeight;
        gph = gph/1000;
		o3_parital_pressure = csv_data.O3PartialPressure;
		p = csv_data.Pressure;
		o3_int = interp1(gph, o3_parital_pressure, altitude_int, 'linear');
		p_int = interp1(gph, p, altitude_int, 'linear');
		o3_vmr = 1000 * (o3_int)./(100 * p_int);

		o3_pc = NaN(size(altitude_int));

		for i = 1:size(altitude_int,2)-1
			o3_pc(i+1) = (10*r*t0)*(0.5*(o3_vmr(i)+o3_vmr(i+1))*(p_int(i)-p_int(i+1)))/(g*p0);
		end

		o3_pc_smt = NaN(size(altitude));

		for i = 1:size(altitude,2)
			o3_pc_smt(i) = nansum(o3_pc(10*(i-1)+1:10*i));
		end

		sonde_smoothed_profile = o3_pc_smt*ak';

        save(['/data01/satellite/IASI/O3/PROFILE/METOPB/ozonesonde_poh/ozonesonde_poh_smoothed' date_str '.mat'], 'sonde_smoothed_profile');
		disp(['Saved data for ' date_str]);

    catch
        disp(['No data for same day ' date_str]);
    end 
    
    toc
end





