clear all;
clear model data params options
rng(888);
addpath(fullfile('..', 'mcmcstat'));

input_dir = fullfile('..', 'datasets');
output_dir = fullfile('..', 'results', 'C1');
m = matfile(fullfile(input_dir, 'data_3.mat'));
data = m.data_3;
HBVdata20241 = data.HBVdata20241;

model.ssfun = @mcmc_ss; 

clear all;
clear model data params options
rng(888);
input_dir = fullfile('..', 'datasets');
output_dir = fullfile('..', 'results', 'C1');
m = matfile(fullfile(input_dir, 'data_3.mat'));
data = m.data_3;
HBVdata20241 = data.HBVdata20241;

model.ssfun = @mcmc_ss; 
params = {
    {'a0', 0.25, 0, 1}
    {'a1', 0.25, 0, 1}
    {'a2', 0.25, 0, 1}
    {'a3', 0.25, 0, 1}
    {'a4', 0.25, 0, 1}
    {'a5', 0.25, 0, 1}
    {'a6', 0.25, 0, 1}
    {'a7', 0.30, 0, 1}
    {'a8', 0.30, 0, 1}
    {'a9', 0.30, 0, 1}
    {'a10', 0.30, 0, 1}
    {'a11', 0.30, 0, 1}
    {'a12', 0.30, 0, 1}
    {'a13', 0.30, 0, 1}
    {'a14', 0.30, 0, 1} 
    {'a15', 0.30, 0, 1}
    {'a16', 0.30, 0, 1}
    {'a17', 0.30, 0, 1}
    {'a18', 0.30, 0, 1}
    {'a19', 0.30, 0, 1}
    {'a20', 0.30, 0, 1}
    {'m1', 0.3, 0, Inf} 
    {'m2', 0.08,  0, Inf}
    {'m3', 55, 0, Inf}
    {'m4', 0.01,0, Inf}
    {'m5', 0.4, 0, Inf}
    {'m6', 0.08,  0, Inf}
    {'m7', 60, 0, Inf}
    {'m8', 0.01, 0, Inf}
    {'u1', 65, 0, Inf}
    {'u2', 30, 0, Inf}
    {'totalcc1992', 2500000, 0, 5000000}
    {'totaldc1992', 800000,  0, 1500000}
    {'totalhcc1992',1200000, 0, 3000000}
    };     

model.S20 = 2;
model.N0 = 10;
model.N = 180;

options.nsimu = 500000; 
options.method = 'dram'; 
options.verbosity = 0;
[results, chain1, s2chain1] = mcmcrun(model, data, params, options);
global zlthbvtmchb zlthbvtmdc zlthbvtmhcc hccdeathmat hccnewmat zlthbvtmcc
options.nsimu = 500000;
[results, chain2, s2chain2] = mcmcrun(model, data, params, options, results);

chain_burnt = [chain2];
[stats] = chainstats(chain_burnt, results);

save(fullfile(output_dir, 'MCMC_results_C1.mat'), ...
    'results', 'chain1', 's2chain1', 'chain2', 's2chain2', 'chain_burnt', 'stats', '-v7.3');
load(fullfile(output_dir, 'MCMC_results_C1.mat'), ...
    'results', 'chain1', 's2chain1', 'chain2', 's2chain2', 'chain_burnt', 'stats');
writematrix(zlthbvtmhcc, fullfile(output_dir, 'zlthbvtmhcc.csv'));
writematrix(hccnewmat, fullfile(output_dir, 'hccnewmat.csv'));
writematrix(hccdeathmat, fullfile(output_dir, 'hccdeathmat.csv'));

figure('Position', [100, 100, 1200, 800]);
mcmcplot(chain_burnt(:, 1:12), [], results, 'chainpanel');

figure('Position', [150, 150, 1200, 800]);
mcmcplot(chain_burnt(:, 13:24), [], results, 'chainpanel');

figure('Position', [200, 200, 1200, 800]);
mcmcplot(chain_burnt(:, 25:end), [], results, 'chainpanel');

modelfun = @(data,th) mcmc_model(th, th(end-2:end), data, 23);
nsample = 5000;
out = mcmcpred(results, chain_burnt, [], data, modelfun, nsample);

years = 1:19;  
obs = HBVdata20241(1:19, 36:46); 
xticklabel = {'<1','1–4','5–9','10–14','15–19','20–24','25–29', ...
              '30–34','35–39','40–44','45–49','50–54','55–59', ...
              '60–64','65–69','70–74','75–79','80–84','>=85'};
figure('Position', [100, 100, 1400, 800]);
for j = 1:11
    pred = out.predlims{1,1}{1:j};  
    subplot(3, 4, j); hold on;
    fill([years, fliplr(years)], ...
         [pred(1,:), fliplr(pred(9,:))], ...
         [1 0.8 0.8], 'EdgeColor', 'none');   

    plot(years, (pred(1,:)+pred(9,:))/2, 'r-', 'LineWidth', 1.5); 
    plot(years, obs(:,j), 'ko', 'MarkerFaceColor', 'k');

    xlim([1, 19]);
    ylim([0, 120]);
    yticks(0:20:110); 

    xticks(1:19);
    xticklabels(xticklabel);
    xtickangle(45);
    set(gca, 'FontSize', 7); 
    grid off; 
end

tmp1 = num2cell(stats(22:29)); [m1, m2, m3, m4, m5, m6, m7, m8] = tmp1{:};
x = (0:100)';
f1 = @(x) m1 ./ (1 + exp(-m2 * (x - m3))) + m4;  
f2 = @(x) m5 ./ (1 + exp(-m6 * (x - m7))) + m8; 
x = (0:100)';y1 = f1(x);y2 = f2(x);
figure;
plot(x, y1, 'b-', 'LineWidth', 2); hold on;
plot(x, y2, 'r-', 'LineWidth', 2); hold on;
ylim([0, 1.0]);
grid on;

median_vals = prctile(chain_burnt, 50)';  
lower_bounds = prctile(chain_burnt, 2.5)';
upper_bounds = prctile(chain_burnt, 97.5)';
std_vals = std(chain_burnt)';
tau = stats(:, 4); 
ess = size(chain_burnt, 1) ./ tau; 
geweke = stats(:, 5); 
ci_table = [median_vals, lower_bounds, upper_bounds, std_vals, ess, geweke];
writematrix(ci_table, fullfile(output_dir, 'ci_table.csv'));

best_theta = median_vals'; 
y0_params = best_theta(32:34);

pred_mortality = mcmc_model(best_theta, y0_params, data, 23);
obs_mortality = HBVdata20241(1:19, 36:46);

residuals = obs_mortality - pred_mortality;

rmse = sqrt(mean(residuals(:).^2));

valid_idx = obs_mortality(:) > 0;
mape = mean(abs(residuals(valid_idx) ./ obs_mortality(valid_idx))) * 100;

SS_res = sum(residuals(:).^2);
SS_tot = sum((obs_mortality(:) - mean(obs_mortality(:))).^2);
R2 = 1 - (SS_res / SS_tot);

param_corr_matrix = corrcoef(chain_burnt);

param_names = cellfun(@(x) x{1}, params, 'UniformOutput', false);

figure('Position', [100, 100, 800, 600]);
imagesc(param_corr_matrix);
colorbar;
colormap(jet); 
caxis([-1 1]); 
xticks(1:length(param_names));
xticklabels(param_names);
xtickangle(45);
yticks(1:length(param_names));
yticklabels(param_names);
set(gca, 'FontSize', 8);

T_corr = array2table(param_corr_matrix, 'VariableNames', param_names, 'RowNames', param_names);
writetable(T_corr, fullfile(output_dir, 'Parameter_Correlation_Matrix.csv'), 'WriteRowNames', true);

figure('Position', [100, 100, 800, 600]);
imagesc(param_corr_matrix);
colorbar;
colormap(jet); 
caxis([-1 1]); 
xticks(1:length(param_names));
xticklabels(param_names);
xtickangle(45);
yticks(1:length(param_names));
yticklabels(param_names);
set(gca, 'FontSize', 8);

writetable(T_corr, fullfile(output_dir, 'Parameter_Correlation_Matrix.csv'), 'WriteRowNames', true);

param_names = cellfun(@(x) x{1}, params, 'UniformOutput', false);
figure('Position', [50, 50, 1600, 1000]);
mcmcplot(chain_burnt, [], results, 'chainpanel');

axes_list = findall(gcf, 'type', 'axes');
for i = 1:length(axes_list)
    ax = axes_list(i);
    if ~isempty(ax.Title.String)
        p_idx = find(strcmp(param_names, ax.Title.String)); 
        
        if ~isempty(p_idx)
            yl = ylim(ax);
           
            y_min = params{p_idx}{3}; 
            y_max = params{p_idx}{4};
            
            if strcmp(param_names{p_idx}, 'totalcc1992')
                y_max = 6000000;
            end
            
            if isinf(y_max)
                y_max = yl(2); 
            end 
            if isinf(y_min)
                y_min = yl(1); 
            end
            
            ylim(ax, [y_min, y_max]); 
            xlim(ax, [1, size(chain_burnt, 1)]); 
            
            set(ax, 'FontSize', 7);
        end
    end
end