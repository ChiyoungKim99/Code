avogadro_number = 6.0221e23;
du_factor = 2.69e16;

file_path = ('/data01/satellite/TROPOMI/O3/TCO/');
mat_file_path = ('/data01/satellite/TROPOMI/O3/TCO/mat/');

start_date = datenum('20200101', 'yyyymmdd');
end_date = datenum('20201231', 'yyyymmdd');
date_range = start_date:end_date;

for t = 1:numel(date_range)
    tic
    date_str = datestr(date_range(t), 'yyyymmdd');
    file_pattern = [file_path 'S5P_RPRO_L2__O3_____' date_str '*'];
    tropomi_files = dir(file_pattern);
    
    tco = [];
    lat = [];
    lon = [];
    
    for file_idx = 1:numel(tropomi_files)
        file_tco = fullfile(file_path, tropomi_files(file_idx).name);
        tco_temp = squeeze(ncread(file_tco, '/PRODUCT/ozone_total_vertical_column')) * avogadro_number / (1e4 * du_factor);
        lat_temp = squeeze(ncread(file_tco, '/PRODUCT/latitude'));
        lon_temp = squeeze(ncread(file_tco, '/PRODUCT/longitude'));
        qa_temp = squeeze(ncread(file_tco, '/PRODUCT/qa_value'));
        tco_temp(qa_temp < 0.5) = [];
        lat_temp(qa_temp < 0.5) = [];
        lon_temp(qa_temp < 0.5) = [];
    end
    % Create a table from the data and save it as a CSV file
    data_table = table(tco_temp(:), lat_temp(:), lon_temp(:), 'VariableNames', {'tco', 'lat', 'lon'});
	csv_file_name = fullfile(mat_file_path, ['TROPOMI_' date_str '.csv']);
	writetable(data_table, csv_file_name);
    toc
end