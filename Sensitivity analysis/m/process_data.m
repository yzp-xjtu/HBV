clear model data params options

global SENS_CFG;
if isempty(SENS_CFG)
    SENS_CFG.dir_C = fullfile('..', 'results', 'C1');
end

input_dir1 = fullfile('..', 'raw_data');
input_dir2 = SENS_CFG.dir_C;
output_dir = SENS_CFG.dir_C; 

file_path1 = fullfile(input_dir1, 'HBVdata20241.xlsx'); 
file_path2 = fullfile(input_dir1, 'HBVdata20242.xlsx'); 
file_path3 = fullfile(input_dir1, 'HBVdata20243.xlsx'); 
file_path5 = fullfile(input_dir2, '1_zlthbvtm4101.csv'); 
file_path4 = fullfile(input_dir2, '1_zlthbvtm2101.csv'); 
file_path6 = fullfile(input_dir2, '2_zlthbvtm2101.csv'); 
file_path7 = fullfile(input_dir2, '3_zlthbvtm2101.csv'); 

HBVdata20241 = readmatrix(file_path1, 'NumHeaderLines', 1);  
HBVdata20242 = readmatrix(file_path2);
HBVdata20243 = readmatrix(file_path3);
HBVdata20245 = readmatrix(file_path5);
HBVdata20244 = readmatrix(file_path4);
HBVdata20246 = readmatrix(file_path6);
HBVdata20247 = readmatrix(file_path7);

save(fullfile(output_dir, 'HBVdata20241.mat'), 'HBVdata20241');
save(fullfile(output_dir, 'HBVdata20242.mat'), 'HBVdata20242');
save(fullfile(output_dir, 'HBVdata20243.mat'), 'HBVdata20243');
save(fullfile(output_dir, 'HBVdata20245.mat'), 'HBVdata20245');

data1 = struct();
data1.HBVdata20241 = HBVdata20241;
data1.HBVdata20242 = HBVdata20242;
data1.HBVdata20243 = HBVdata20243;
save(fullfile(output_dir, 'data1.mat'), 'data1');

data_1 = struct();
data_1.HBVdata20241 = HBVdata20241;
data_1.HBVdata20242 = HBVdata20242;
data_1.HBVdata20243 = HBVdata20243;
data_1.HBVdata20244 = HBVdata20244;
data_1.HBVdata20245 = HBVdata20245;
save(fullfile(output_dir, 'data_1.mat'), 'data_1');

data_2 = struct();
data_2.HBVdata20241 = HBVdata20241;
data_2.HBVdata20242 = HBVdata20242;
data_2.HBVdata20243 = HBVdata20243;
data_2.HBVdata20244 = HBVdata20246;
data_2.HBVdata20245 = HBVdata20245;
save(fullfile(output_dir, 'data_2.mat'), 'data_2');

data_3 = struct();
data_3.HBVdata20241 = HBVdata20241;
data_3.HBVdata20242 = HBVdata20242;
data_3.HBVdata20243 = HBVdata20243;
data_3.HBVdata20244 = HBVdata20247;
data_3.HBVdata20245 = HBVdata20245;
save(fullfile(output_dir, 'data_3.mat'), 'data_3');