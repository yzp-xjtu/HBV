clear all; clc;
addpath(fullfile('..', 'sensitivity_analysis_code'));
global SENS_CFG;

base_sens_dir = fullfile('..', 'sensitivity_analysis_results');
if ~exist(base_sens_dir, 'dir'), mkdir(base_sens_dir); end

run_scenario('Vaccine_Protection', {'0.75', '0.8', '0.9', '0.95'}, 'full');

run_scenario('Clearance_Rate', {'0.005', '0.008', '0.012', '0.015'}, 'full');

run_scenario('Transmission', {'lower', 'upper'}, 'full');

run_scenario('Intranterine', {'0.7', '0.75', '0.8', '0.85'}, 'full');

v_names = {'v1', 'v2', 'v3', 'v4', 'v5'};
v_mults = [0.8, 0.9, 1.1, 1.2];
for v_idx = 1:5
    v_str = v_names{v_idx}; 
    for m_idx = 1:length(v_mults)
        val_str = sprintf('%.1f', v_mults(m_idx));
        run_scenario(sprintf('Prog_Rate_%s', v_str), {val_str}, 'hbv_only', v_str, v_mults(m_idx));
    end
end

run_scenario('Prog_Rate_Joint', {'random_20'}, 'hbv_only_joint');

run_scenario('Discount_Rate', {'0.00', '0.03', '0.05', '0.07'}, 'economic_only');

disp('==================================================');
disp('All sensitivity analysis tasks completed!');

function run_scenario(param_name, values, run_type, varargin)
    global SENS_CFG;
    base_original_C = fullfile('..', 'results', 'C1');
    base_original_HBV = fullfile('..', 'results', 'HBV');
    base_sens_dir = fullfile('..', 'sensitivity_analysis_results');

    for v = 1:length(values)
        val = values{v};
        
        out_root = fullfile(base_sens_dir, sprintf('%s(%s)', param_name, val));
        out_C = fullfile(out_root, 'C1');
        out_HBV = fullfile(out_root, 'HBV');
        
        if ~exist(out_C, 'dir'), mkdir(out_C); end
        if ~exist(out_HBV, 'dir'), mkdir(out_HBV); end
        
        SENS_CFG.p_vaccine = 0.85; 
        SENS_CFG.clearance_rate = 'random';
        SENS_CFG.transm_bound = 'random'; 
        SENS_CFG.intranterine = 0.8;
        SENS_CFG.v_target = 'none'; 
        SENS_CFG.v_val = 1.0;
        SENS_CFG.discount_rate = 0.03;
        SENS_CFG.dir_C = out_C;
        SENS_CFG.dir_HBV = out_HBV;
        
        switch param_name
            case 'Vaccine_Protection'
                SENS_CFG.p_vaccine = str2double(val);
            case 'Clearance_Rate'
                SENS_CFG.clearance_rate = val;
            case 'Transmission'
                SENS_CFG.transm_bound = val;
            case 'Intranterine'
                SENS_CFG.intranterine = str2double(val);
            case 'Discount_Rate'
                SENS_CFG.discount_rate = str2double(val);
            case 'Prog_Rate_Joint'
                SENS_CFG.v_target = 'joint';
            otherwise
                if contains(param_name, 'Prog_Rate_')
                    SENS_CFG.v_target = varargin{1}; 
                    SENS_CFG.v_val = varargin{2};    
                end
        end
        
        c1_file_check  = fullfile(out_C, '3_zlthbvtm2101.csv'); 
        hbv_file_check = fullfile(out_HBV, '3_zlthbvtmhcc103.csv');
        d_file_check   = fullfile(out_C, '7_2024.xlsx');
        
        c1_done  = exist(c1_file_check, 'file');
        hbv_done = exist(hbv_file_check, 'file');
        d_done   = exist(d_file_check, 'file');

        switch run_type
            case 'full'
                if d_done && hbv_done && c1_done
                    cleanup_large_mats(out_C, out_HBV); 
                    continue;
                end
                
                if ~c1_done
                    C1_2; C2_2; C3_3; process_data_2;
                end
                
                if ~hbv_done
                    HBV1_2; HBV2_2; HBV3_2;
                end
                
                if ~d_done
                    d;
                end
                
            case {'hbv_only', 'hbv_only_joint'}
                if d_done && hbv_done
                    cleanup_large_mats(out_C, out_HBV);
                    continue;
                end
                
                copyfile(fullfile(base_original_C, '*.*'), out_C); 
                
                if ~hbv_done
                    HBV1_2; HBV2_2; HBV3_2;
                end
                
                if ~d_done
                    d;
                end
                
            case 'economic_only'
                if d_done
                    cleanup_large_mats(out_C, out_HBV);
                    continue;
                end
                
                copyfile(fullfile(base_original_C, '*.*'), out_C);
                copyfile(fullfile(base_original_HBV, '*.*'), out_HBV);
                d;
        end
        
        if exist(d_file_check, 'file')
            cleanup_large_mats(out_C, out_HBV);
        end
    end
end

function cleanup_large_mats(out_C, out_HBV)
    files_to_delete = {
        fullfile(out_C, 'C1.mat'),
        fullfile(out_C, 'C2.mat'),
        fullfile(out_C, 'C3.mat'),
        fullfile(out_HBV, 'HBV1.mat'),
        fullfile(out_HBV, 'HBV2.mat'),
        fullfile(out_HBV, 'HBV3.mat')
    };
    
    for i = 1:length(files_to_delete)
        if exist(files_to_delete{i}, 'file')
            delete(files_to_delete{i});
        end
    end
end