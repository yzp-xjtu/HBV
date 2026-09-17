%% 数据
clear all;
output_dir_C = fullfile('..', 'results', 'C1');
output_dir_HBV = fullfile('..', 'results', 'HBV');

disp('Loading C scenario data...');
data_C1 = load(fullfile(output_dir_C, 'C1.mat'));
data_C2 = load(fullfile(output_dir_C, 'C2.mat'));
data_C3 = load(fullfile(output_dir_C, 'C3.mat'));

disp('Loading HBV fine-stage data...');
data_HBV1 = load(fullfile(output_dir_HBV, 'HBV1.mat'));
data_HBV2 = load(fullfile(output_dir_HBV, 'HBV2.mat'));
data_HBV3 = load(fullfile(output_dir_HBV, 'HBV3.mat'));

C_total = data_C3.zlthbvtm2(:,1:109);

HBV_total = data_HBV3.zlthbvtmchb(:,1:109,:) + ...
            data_HBV3.zlthbvtmcc(:,1:109,:) + ...
            data_HBV3.zlthbvtmdc(:,1:109,:) + ...
            data_HBV3.zlthbvtmhcc(:,1:109,:);

HBV_total_median = median(HBV_total, 3);
C_total_sum = sum(C_total, 1);
HBV_total_sum = sum(HBV_total_median, 1);
C_greater_HBV = C_total_sum > HBV_total_sum;
[idx_age, idx_year] = find(C_total <= HBV_total_median);

years_1 = 1993:2024;
years_2 = 2025:2100;
years_3 = 1993:2100;

get_yr_idx_C = @(yrs) yrs - 1992 + 1;
get_yr_idx_HBV = @(yrs) yrs - 1993 + 1;

age_group_labels = {"1-4", "5-9", "10-14", "15-19", ...
                    "20-24", "25-29", "30-34", "35-39", ...
                    "40-49", "50-59", "60-69", "70-79", ...
                    "80-89", "90-100"};

age_groups_idx = {
    2:5,      
    6:10,     
    11:15,    
    16:20,    
    21:25,    
    26:30,    
    31:35,    
    36:40,    
    41:50,    
    51:60,    
    61:70,    
    71:80,    
    81:90,    
    91:100    
};

disp('Calculating Table 1...');
survey_years = [2006, 2014, 2020];
t1_rows = {};
for y = survey_years
    y_idx = get_yr_idx_C(y);
    for g = 1:length(age_group_labels)
        ages = age_groups_idx{g};
        r1 = squeeze(sum(data_C1.zlthbvtm2(ages+1, y_idx, :), 1) ./ sum(data_C1.zlthbvtm4(ages+1, y_idx, :), 1)) * 100;
        r2 = squeeze(sum(data_C2.zlthbvtm2(ages+1, y_idx, :), 1) ./ sum(data_C2.zlthbvtm4(ages+1, y_idx, :), 1)) * 100;
        r3 = squeeze(sum(data_C3.zlthbvtm2(ages+1, y_idx, :), 1) ./ sum(data_C3.zlthbvtm4(ages+1, y_idx, :), 1)) * 100;
        
        fmt_rate = @(v) sprintf('%.2f (%.2f, %.2f)', prctile(v, 50), prctile(v, 2.5), prctile(v, 97.5));
        pct1 = (r1 - r2) ./ r1 * 100;
        pct2 = (r1 - r2) ./ (r1 - r3) * 100;
        fmt_pct = @(v) sprintf('%.2f (%.2f, %.2f)', prctile(v, 50), prctile(v, 2.5), prctile(v, 97.5));
        
        t1_rows(end+1, :) = {y, char(age_group_labels{g}), fmt_rate(r1), fmt_rate(r2), fmt_rate(r3), fmt_pct(pct1), fmt_pct(pct2)};
    end
end
T1 = cell2table(t1_rows, 'VariableNames', {'year', 'ages', 'HBV1', 'HBV2', 'HBV3', 'Reduced_by_newborn_vaccination_pct', 'Effect_within_comprehensive_interventions_pct'});
writetable(T1, fullfile(output_dir_C, '1.xlsx'));

disp('Calculating Table 2...');
t2_rows = {};
for y = survey_years
    y_idx = get_yr_idx_C(y);
    ages = 1:59;
    r1 = squeeze(sum(data_C1.zlthbvtm2(ages+1, y_idx, :), 1) ./ sum(data_C1.zlthbvtm4(ages+1, y_idx, :), 1)) * 100;
    r2 = squeeze(sum(data_C2.zlthbvtm2(ages+1, y_idx, :), 1) ./ sum(data_C2.zlthbvtm4(ages+1, y_idx, :), 1)) * 100;
    r3 = squeeze(sum(data_C3.zlthbvtm2(ages+1, y_idx, :), 1) ./ sum(data_C3.zlthbvtm4(ages+1, y_idx, :), 1)) * 100;
    
    fmt_rate = @(v) sprintf('%.2f (%.2f, %.2f)', prctile(v, 50), prctile(v, 2.5), prctile(v, 97.5));
    pct1 = (r1 - r2) ./ r1 * 100;
    pct2 = (r1 - r2) ./ (r1 - r3) * 100;
    fmt_pct = @(v) sprintf('%.2f (%.2f, %.2f)', prctile(v, 50), prctile(v, 2.5), prctile(v, 97.5));
    
    t2_rows(end+1, :) = {y, fmt_rate(r1), fmt_rate(r2), fmt_rate(r3), fmt_pct(pct1), fmt_pct(pct2)};
end
T2 = cell2table(t2_rows, 'VariableNames', {'year', 'HBV1', 'HBV2', 'HBV3', 'Reduced_by_newborn_vaccination_pct', 'Effect_within_comprehensive_interventions_pct'});
writetable(T2, fullfile(output_dir_C, '2.xlsx'));

disp('Calculating Table 3...');
p_idx_list = {get_yr_idx_HBV(years_1), get_yr_idx_HBV(years_2), get_yr_idx_HBV(years_3)};
periods = {'period1_1993_2024', 'period2_2025_2100', 'period3_1993_2100'};
target_vars_t3 = {
    'C',   'zlthbvtm5',   'C1';
    'CHB', 'zlthbvtm6',   'C1';
    'CC',  'ccnewmat',    'HBV';
    'DC',  'dcnewmat',    'HBV';
    'HCC', 'hccnewmat',   'HBV'
};
t3_rows = {};
for p = 1:length(periods)
    yr_idx = p_idx_list{p};
    row_data = {periods{p}};
    
    for v = 1:size(target_vars_t3, 1)
        v_name = target_vars_t3{v, 2};
        d_type = target_vars_t3{v, 3};
        
        s1 = zeros(5000, 1); s2 = zeros(5000, 1); s3 = zeros(5000, 1);
        for y = yr_idx
            if strcmp(d_type, 'C1')
                s1 = s1 + squeeze(sum(data_C1.(v_name)(1:101, y, :), 1));
                s2 = s2 + squeeze(sum(data_C2.(v_name)(1:101, y, :), 1));
                s3 = s3 + squeeze(sum(data_C3.(v_name)(1:101, y, :), 1));
            else
                s1 = s1 + squeeze(sum(data_HBV1.(v_name)(1:101, y, :), 1));
                s2 = s2 + squeeze(sum(data_HBV2.(v_name)(1:101, y, :), 1));
                s3 = s3 + squeeze(sum(data_HBV3.(v_name)(1:101, y, :), 1));
            end
        end
        
        red_vac = s1 - s2;
        red_comp = s1 - s3;
        vac_pct = (s1 - s2) ./ s1 * 100;
        vac_eff = (s1 - s2) ./ (s1 - s3) * 100;
        
        fmt_val = @(v_vec) sprintf('%.2f (%.2f, %.2f)', prctile(v_vec, 50)/1e6, prctile(v_vec, 2.5)/1e6, prctile(v_vec, 97.5)/1e6);
        fmt_pct = @(v_vec) sprintf('%.2f (%.2f, %.2f)', prctile(v_vec, 50), prctile(v_vec, 2.5), prctile(v_vec, 97.5));
        
        row_data = [row_data, {fmt_val(s1), fmt_val(s2), fmt_val(s3), ...
                              fmt_val(red_vac), fmt_val(red_comp), ...
                              fmt_pct(vac_pct), fmt_pct(vac_eff)}];
    end
    t3_rows(end+1, :) = row_data;
end
var_names_t3 = {'period'};
groups_t3 = {'C', 'CHB', 'CC', 'DC', 'HCC'};
for v = 1:length(groups_t3)
    nm = groups_t3{v};
    var_names_t3 = [var_names_t3, {[nm '_1'], [nm '_2'], [nm '_3'], ...
                                  [nm '_Vaccine_Reduction'], [nm '_Comprehensive_Reduction'], ...
                                  [nm '_Vaccine_Reduction_%'], [nm '_Vaccine_Effect']}];
end
T3 = cell2table(t3_rows, 'VariableNames', var_names_t3);
writetable(T3, fullfile(output_dir_C, '3.xlsx'));

disp('Calculating Table 4...');
target_vars_t4 = {
    'hbvdeath', 'HBVdeath_all', 'C1';
    'dcdeath',  'dcdeathmat',   'HBV';
    'hccdeath', 'hccdeathmat',  'HBV'
};
t4_rows = {};
for p = 1:length(periods)
    yr_idx = p_idx_list{p};
    row_data = {periods{p}};
    
    for v = 1:size(target_vars_t4, 1)
        v_name = target_vars_t4{v, 2};
        d_type = target_vars_t4{v, 3};
        
        s1 = zeros(5000, 1); s2 = zeros(5000, 1); s3 = zeros(5000, 1);
        for y = yr_idx
            if strcmp(d_type, 'C1')
                s1 = s1 + squeeze(sum(data_C1.(v_name)(1:101, y, :), 1));
                s2 = s2 + squeeze(sum(data_C2.(v_name)(1:101, y, :), 1));
                s3 = s3 + squeeze(sum(data_C3.(v_name)(1:101, y, :), 1));
            else
                s1 = s1 + squeeze(sum(data_HBV1.(v_name)(1:101, y, :), 1));
                s2 = s2 + squeeze(sum(data_HBV2.(v_name)(1:101, y, :), 1));
                s3 = s3 + squeeze(sum(data_HBV3.(v_name)(1:101, y, :), 1));
            end
        end
        
        red_vac = s1 - s2;
        red_comp = s1 - s3;
        vac_pct = (s1 - s2) ./ s1 * 100;
        vac_eff = (s1 - s2) ./ (s1 - s3) * 100;
        
        fmt_val = @(v_vec) sprintf('%.2f (%.2f, %.2f)', prctile(v_vec, 50)/1e6, prctile(v_vec, 2.5)/1e6, prctile(v_vec, 97.5)/1e6);
        fmt_pct = @(v_vec) sprintf('%.2f (%.2f, %.2f)', prctile(v_vec, 50), prctile(v_vec, 2.5), prctile(v_vec, 97.5));
        
        row_data = [row_data, {fmt_val(s1), fmt_val(s2), fmt_val(s3), ...
                              fmt_val(red_vac), fmt_val(red_comp), ...
                              fmt_pct(vac_pct), fmt_pct(vac_eff)}];
    end
    t4_rows(end+1, :) = row_data;
end
var_names_t4 = {'Period'};
groups_t4 = {'hbvdeath', 'dcdeath', 'hccdeath'};
for v = 1:length(groups_t4)
    nm = groups_t4{v};
    var_names_t4 = [var_names_t4, {[nm '_1'], [nm '_2'], [nm '_3'], ...
                                  [nm '_Vaccine_Reduction'], [nm '_Comprehensive_Reduction'], ...
                                  [nm '_Vaccine_Reduction_%'], [nm '_Vaccine_Effect']}];
end
T4 = cell2table(t4_rows, 'VariableNames', var_names_t4);
writetable(T4, fullfile(output_dir_C, '4.xlsx'));

disp('Calculating Table 5...');
years_t5 = [2006, 2014, 2020, 2030, 2050, 2100];
t5_rows = {};

for y = years_t5
    y_idx_C = get_yr_idx_C(y);
    y_idx_HBV = get_yr_idx_HBV(y);
    
    for g = 1:length(age_group_labels)
        ages = age_groups_idx{g};
        
        metrics = {'cnew', 'Rate', 'HBVdeath'};
        for m_idx = 1:length(metrics)
            met = metrics{m_idx};
            
            if strcmp(met, 'cnew')
                v1_vec = squeeze(sum(data_C1.zlthbvtm5(ages+1, y_idx_HBV, :), 1));
                v2_vec = squeeze(sum(data_C2.zlthbvtm5(ages+1, y_idx_HBV, :), 1));
                v3_vec = squeeze(sum(data_C3.zlthbvtm5(ages+1, y_idx_HBV, :), 1));
            elseif strcmp(met, 'Rate')
                v1_vec = squeeze(sum(data_C1.zlthbvtm2(ages+1, y_idx_C, :), 1) ./ sum(data_C1.zlthbvtm4(ages+1, y_idx_C, :), 1)) * 100;
                v2_vec = squeeze(sum(data_C2.zlthbvtm2(ages+1, y_idx_C, :), 1) ./ sum(data_C2.zlthbvtm4(ages+1, y_idx_C, :), 1)) * 100;
                v3_vec = squeeze(sum(data_C3.zlthbvtm2(ages+1, y_idx_C, :), 1) ./ sum(data_C3.zlthbvtm4(ages+1, y_idx_C, :), 1)) * 100;
            elseif strcmp(met, 'HBVdeath')
                v1_vec = squeeze(sum(data_HBV1.dcdeathmat(ages+1, y_idx_HBV, :) + data_HBV1.hccdeathmat(ages+1, y_idx_HBV, :), 1));
                v2_vec = squeeze(sum(data_HBV2.dcdeathmat(ages+1, y_idx_HBV, :) + data_HBV2.hccdeathmat(ages+1, y_idx_HBV, :), 1));
                v3_vec = squeeze(sum(data_HBV3.dcdeathmat(ages+1, y_idx_HBV, :) + data_HBV3.hccdeathmat(ages+1, y_idx_HBV, :), 1));
            end
            
            diff1_vec = v1_vec - v2_vec;
            diff2_vec = v1_vec - v3_vec;
            p1_vec = (v1_vec - v2_vec) ./ v1_vec * 100;
            p2_vec = (v1_vec - v2_vec) ./ (v1_vec - v3_vec) * 100;

            p1_valid = p1_vec(isfinite(p1_vec));
            p2_valid = p2_vec(isfinite(p2_vec));

            p1_str = sprintf('%.2f (%.2f, %.2f)', ...
                prctile(p1_valid, 50), ...
                prctile(p1_valid, 2.5), ...
                prctile(p1_valid, 97.5));

            p2_str = sprintf('%.2f (%.2f, %.2f)', ...
                prctile(p2_valid, 50), ...
                prctile(p2_valid, 2.5), ...
                prctile(p2_valid, 97.5));
            
            fmt_val = @(v_vec) sprintf('%.2f (%.2f, %.2f)', prctile(v_vec, 50), prctile(v_vec, 2.5), prctile(v_vec, 97.5));
            t5_rows(end+1, :) = {char(age_group_labels{g}), met, y, ... 
                     fmt_val(v1_vec), fmt_val(v2_vec), fmt_val(v3_vec), ... 
                     fmt_val(diff1_vec), fmt_val(diff2_vec), ... 
                     p1_str, p2_str};
        end
    end
end

T5 = cell2table(t5_rows, 'VariableNames', {'Age', 'Metric', 'Year', 'HBV1', 'HBV2', 'HBV3', 'Diff1', 'Diff2', 'Percent1', 'Percent2'});
writetable(T5, fullfile(output_dir_C, '5.xlsx'));

t6_rows = {};
groups_t6 = {
    'C',        'zlthbvtm5',   'C1';
    'CHB',      'zlthbvtm6',   'C1';
    'CC',       'ccnewmat',    'HBV';
    'DC',       'dcnewmat',    'HBV';
    'HCC',      'hccnewmat',   'HBV';
    'dcdeath',  'dcdeathmat',  'HBV';
    'hccdeath', 'hccdeathmat', 'HBV'
};
for p = 1:length(periods)
    yr_idx = p_idx_list{p};
    
    for g = 1:length(age_group_labels)
        ages = age_groups_idx{g};
        row_data = {periods{p}, char(age_group_labels{g})};
        
        for v = 1:size(groups_t6, 1)
            v_name = groups_t6{v, 2};
            d_type = groups_t6{v, 3};
            
            s1 = zeros(5000, 1); s2 = zeros(5000, 1); s3 = zeros(5000, 1);
            for y = yr_idx
                if strcmp(d_type, 'C1')
                    s1 = s1 + squeeze(sum(data_C1.(v_name)(ages+1, y, :), 1));
                    s2 = s2 + squeeze(sum(data_C2.(v_name)(ages+1, y, :), 1));
                    s3 = s3 + squeeze(sum(data_C3.(v_name)(ages+1, y, :), 1));
                else
                    s1 = s1 + squeeze(sum(data_HBV1.(v_name)(ages+1, y, :), 1));
                    s2 = s2 + squeeze(sum(data_HBV2.(v_name)(ages+1, y, :), 1));
                    s3 = s3 + squeeze(sum(data_HBV3.(v_name)(ages+1, y, :), 1));
                end
            end
            
            vac_pct = (s1 - s2) ./ s1 * 100;
            vac_eff = (s1 - s2) ./ (s1 - s3) * 100;
            
            fmt_pct = @(v_vec) sprintf('%.2f (%.2f, %.2f)', prctile(v_vec, 50), prctile(v_vec, 2.5), prctile(v_vec, 97.5));
            row_data = [row_data, {fmt_pct(vac_pct), fmt_pct(vac_eff)}];
        end
        t6_rows(end+1, :) = row_data;
    end
end

var_names_t6 = {'Period', 'age'};
for v = 1:size(groups_t6, 1)
    nm = groups_t6{v, 1};
    var_names_t6 = [var_names_t6, {[nm '_Vaccine_Reduction_%'], [nm '_Vaccine_Contribution_Rate']}];
end

T6 = cell2table(t6_rows, 'VariableNames', var_names_t6);
writetable(T6, fullfile(output_dir_C, '6.xlsx'));
disp('Table 6 calculated!');

disp('Calculating Table 7...');
exchange_rate = 6.73;  
discount_rate = 0.03;  
base_year = 1992;      

cost_dir_CHB = 15883.53 + 2980.93; 
cost_dir_CC  = 25540.83;
cost_dir_DC  = 33368.82;
cost_dir_HCC = 45058.51;

cost_ind_CHB = 4592.36;
cost_ind_CC  = 7374.05;
cost_ind_DC  = 10543.60;
cost_ind_HCC = 17452.10;

prod_loss_death = [repmat(1822577.87, 1, 21), repmat(1311194.73, 1, 20), repmat(355737.10, 1, 60)]';

dw_CHB = 0.000;  
dw_CC  = 0.000;  
dw_DC  = 0.180;  
dw_HCC = 0.450;  

SLE = [88.8718951, repmat(88.00051053, 1, 4), repmat(84.03008056, 1, 5), ...
       repmat(79.04633476, 1, 5), repmat(74.0665492, 1, 5), repmat(69.10756792, 1, 5), ...
       repmat(64.14930031, 1, 5), repmat(59.1962771, 1, 5), repmat(54.25261364, 1, 5), ...
       repmat(49.31739311, 1, 5), repmat(44.43332057, 1, 5), repmat(39.63473787, 1, 5), ...
       repmat(34.91488095, 1, 5), repmat(30.25343822, 1, 5), repmat(25.68089534, 1, 5), ...
       repmat(21.28820012, 1, 5), repmat(17.10351469, 1, 5), repmat(13.23872477, 1, 5), ...
       repmat(9.990181244, 1, 5), repmat(7.617724915, 1, 5), repmat(5.922359078, 1, 6)]';

try
    vprice = readmatrix(fullfile('..', 'datasets', 'vaccine.csv'));
    vaccine_years = 1993 : (1993 + size(vprice, 1) - 1);
    
    admin_fee = 23.30;
    other_fee = 41.96 * 1.5; 
    
    proc_cost_sum = 0; admin_cost_sum = 0; other_cost_sum = 0;
    for vy = 1:length(vaccine_years)
        cur_y = vaccine_years(vy);
        df_v = 1 / ((1 + discount_rate)^(cur_y - base_year)); 
        
        cohort   = vprice(vy, 9);
        coverage = vprice(vy, 8);
        price    = vprice(vy, 5);
        
        doses = cohort * coverage * 3;
        proc_cost_sum  = proc_cost_sum  + (doses * price) * df_v;
        admin_cost_sum = admin_cost_sum + (doses * admin_fee) * df_v;
        other_cost_sum = other_cost_sum + (doses * other_fee) * df_v;
    end
    
    Vac_Proc_Billion  = proc_cost_sum  / (exchange_rate * 1e9);
    Vac_Admin_Billion = admin_cost_sum / (exchange_rate * 1e9);
    Vac_Other_Billion = other_cost_sum / (exchange_rate * 1e9);
    Vac_Total_Billion = Vac_Proc_Billion + Vac_Admin_Billion + Vac_Other_Billion;
catch
    Vac_Proc_Billion = 1.45; Vac_Admin_Billion = 3.44; Vac_Other_Billion = 9.29;
    Vac_Total_Billion = 14.18;
end

eval_years = 1993:2024; 

acc_Dir_NV = zeros(1, 5000); acc_Dir_V = zeros(1, 5000);
acc_Ind_NV = zeros(1, 5000); acc_Ind_V = zeros(1, 5000);
acc_Prd_NV = zeros(1, 5000); acc_Prd_V = zeros(1, 5000);
acc_DALY_NV = zeros(1, 5000); acc_DALY_V = zeros(1, 5000);

for y = eval_years
    df = 1 / ((1 + discount_rate)^(y - base_year)); 
    
    idx_stock = y - 1992 + 1; 
    idx_flow  = y - 1993 + 1; 
    
    p_CHB_nv = squeeze(data_HBV1.zlthbvtmchb(:, idx_stock, :));
    p_CC_nv  = squeeze(data_HBV1.zlthbvtmcc(:, idx_stock, :));
    p_DC_nv  = squeeze(data_HBV1.zlthbvtmdc(:, idx_stock, :));
    p_HCC_nv = squeeze(data_HBV1.zlthbvtmhcc(:, idx_stock, :));
    
    p_CHB_v  = squeeze(data_HBV2.zlthbvtmchb(:, idx_stock, :));
    p_CC_v   = squeeze(data_HBV2.zlthbvtmcc(:, idx_stock, :));
    p_DC_v   = squeeze(data_HBV2.zlthbvtmdc(:, idx_stock, :));
    p_HCC_v  = squeeze(data_HBV2.zlthbvtmhcc(:, idx_stock, :));
    
    d_DC_nv  = squeeze(data_HBV1.dcdeathmat(:, idx_flow, :));
    d_HCC_nv = squeeze(data_HBV1.hccdeathmat(:, idx_flow, :));
    d_DC_v   = squeeze(data_HBV2.dcdeathmat(:, idx_flow, :));
    d_HCC_v  = squeeze(data_HBV2.hccdeathmat(:, idx_flow, :));
    
    dir_nv = sum(p_CHB_nv)*cost_dir_CHB + sum(p_CC_nv)*cost_dir_CC + sum(p_DC_nv)*cost_dir_DC + sum(p_HCC_nv)*cost_dir_HCC;
    dir_v  = sum(p_CHB_v)*cost_dir_CHB  + sum(p_CC_v)*cost_dir_CC  + sum(p_DC_v)*cost_dir_DC  + sum(p_HCC_v)*cost_dir_HCC;
    
    ind_nv = sum(p_CHB_nv)*cost_ind_CHB + sum(p_CC_nv)*cost_ind_CC + sum(p_DC_nv)*cost_ind_DC + sum(p_HCC_nv)*cost_ind_HCC;
    ind_v  = sum(p_CHB_v)*cost_ind_CHB  + sum(p_CC_v)*cost_ind_CC  + sum(p_DC_v)*cost_ind_DC  + sum(p_HCC_v)*cost_ind_HCC;
    
    prd_nv = prod_loss_death' * (d_DC_nv + d_HCC_nv);
    prd_v  = prod_loss_death' * (d_DC_v + d_HCC_v);
    
    acc_Dir_NV = acc_Dir_NV + dir_nv * df;   acc_Dir_V = acc_Dir_V + dir_v * df;
    acc_Ind_NV = acc_Ind_NV + ind_nv * df;   acc_Ind_V = acc_Ind_V + ind_v * df;
    acc_Prd_NV = acc_Prd_NV + prd_nv * df;   acc_Prd_V = acc_Prd_V + prd_v * df;
    
    yld_nv = sum(p_CHB_nv)*dw_CHB + sum(p_CC_nv)*dw_CC + sum(p_DC_nv)*dw_DC + sum(p_HCC_nv)*dw_HCC;
    yld_v  = sum(p_CHB_v)*dw_CHB  + sum(p_CC_v)*dw_CC  + sum(p_DC_v)*dw_DC  + sum(p_HCC_v)*dw_HCC;
    
    yll_nv = SLE' * (d_DC_nv + d_HCC_nv);
    yll_v  = SLE' * (d_DC_v + d_HCC_v);
    
    acc_DALY_NV = acc_DALY_NV + (yld_nv + yll_nv) * df;
    acc_DALY_V  = acc_DALY_V  + (yld_v  + yll_v) * df;
end

B_factor = exchange_rate * 1e9;
Tot_Cost_NV = (acc_Dir_NV + acc_Ind_NV + acc_Prd_NV) / B_factor;
Tot_Cost_V  = (acc_Dir_V  + acc_Ind_V  + acc_Prd_V)  / B_factor + Vac_Total_Billion;

Tot_DALY_NV = acc_DALY_NV / 1e6;
Tot_DALY_V  = acc_DALY_V  / 1e6;

Averted_Dir = (acc_Dir_NV - acc_Dir_V) / B_factor;
Averted_Ind = (acc_Ind_NV - acc_Ind_V) / B_factor;
Averted_Prd = (acc_Prd_NV - acc_Prd_V) / B_factor;
Averted_TotalCost = Tot_Cost_NV - Tot_Cost_V;
Averted_DALY = Tot_DALY_NV - Tot_DALY_V;

BCR_array = (Tot_Cost_NV - ((acc_Dir_V + acc_Ind_V + acc_Prd_V)/B_factor)) ./ Vac_Total_Billion;
ICER_array = ((Tot_Cost_V - Tot_Cost_NV) * 1e6) ./ ((Tot_DALY_NV - Tot_DALY_V) * 1e6);

fmt = @(v) sprintf('%.2f (%.2f, %.2f)', prctile(v, 50), prctile(v, 2.5), prctile(v, 97.5));

t7_data = {
    'Total costs of vaccination (billion USD)', '', sprintf('%.2f', Vac_Total_Billion), '';
    '  - Vaccine procurement cost', '', sprintf('%.2f', Vac_Proc_Billion), '';
    '  - Administrative service cost', '', sprintf('%.2f', Vac_Admin_Billion), '';
    '  - Other costs', '', sprintf('%.2f', Vac_Other_Billion), '';
    'Total Benefits (billion USD)', fmt(Tot_Cost_NV), fmt(Tot_Cost_V - Vac_Total_Billion), fmt(Averted_TotalCost);
    '  - Direct medical', fmt(acc_Dir_NV/B_factor), fmt(acc_Dir_V/B_factor), fmt(Averted_Dir);
    '  - Indirect costs', fmt(acc_Ind_NV/B_factor), fmt(acc_Ind_V/B_factor), fmt(Averted_Ind);
    '  - Productivity loss', fmt(acc_Prd_NV/B_factor), fmt(acc_Prd_V/B_factor), fmt(Averted_Prd);
    'DALYs (millions)', fmt(Tot_DALY_NV), fmt(Tot_DALY_V), fmt(Averted_DALY);
    'Benefit-Cost Ratio (BCR)', '', '', fmt(BCR_array);
    'ICER (USD per DALY averted)', '', '', fmt(ICER_array)
};

T7 = cell2table(t7_data, 'VariableNames', {'Item', 'Without_vaccination', 'With_vaccination', 'Averted_or_Ratio'});
writetable(T7, fullfile(output_dir_C, sprintf('7_%d.xlsx', eval_years(end))));