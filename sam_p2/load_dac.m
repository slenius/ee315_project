clear all
close all

% load variants
base_dir = '../sam_cadence/simulation/dac_tb/spectre/schematic/';
var_dirs = dir(base_dir);
data_dirs = cell(1);
wu = [];
j = 1;
for i = 3:length(var_dirs)
    n = var_dirs(i).name;
    if strncmpi(n, 'wu=', 2)
      fprintf('%s\n', var_dirs(i).name);
      d = strcat(base_dir, n, '/psf/');
      data_dirs{j} = d;
      wu(j) = str2num(n(4:end));
      j = j + 1;
    end
end


data = cell(j-1,1);

for i = 1:length(data)
    data{i} = struct();
    data{i}.dacout_2_p = cds_srr(data_dirs{i}, 'tran-tran', 'dacout_2_p');
    data{i}.dacout_2_m = cds_srr(data_dirs{i}, 'tran-tran', 'dacout_2_m');
    data{i}.wu = wu(i);
end

save('dacdata_ab_2.mat');