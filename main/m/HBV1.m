clear all;
clear model data params options
rng(888);

input_dir = fullfile('..', 'datasets');
output_dir1 = fullfile('..', 'results', 'C1');
output_dir2 = fullfile('..', 'results', 'HBV');

m = matfile(fullfile(input_dir, 'data_1.mat'));
data = m.data_1;
stats_file = fullfile(output_dir1, 'ci_table.csv'); 
stats = readmatrix(stats_file);

theta_vals = stats(1:34, 1); 
theta_lower = stats(1:31, 2); 
theta_upper = stats(1:31, 3); 
y0_lower = stats(32:34, 2);
y0_upper = stats(32:34, 3);

N = 132;  
num_samples = 5000;  
theta_samples = zeros(31, num_samples);
y0_samples = zeros(3, num_samples);

for i = 1:num_samples
    theta_samples(:,i) = theta_lower + (theta_upper - theta_lower) .* rand(31,1);
    y0_samples(:,i)    = y0_lower + (y0_upper - y0_lower) .* rand(3,1);
end

for m = 1:num_samples
    theta = theta_samples(:,m);
    y0 = y0_samples(:,m);
    HBVdata20241 = data.HBVdata20241;
    HBVdata20242 = data.HBVdata20242;
    HBVdata20243 = data.HBVdata20243;
    HBVdata20244 = data.HBVdata20244;
    HBVdata20245 = data.HBVdata20245;

    tmp1 = num2cell(theta(22:29)); [m1, m2, m3, m4, m5, m6, m7, m8] = tmp1{:};
    tmp2 = num2cell(theta(1:21)); [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20] = tmp2{:};
    tmp3 = num2cell(theta_vals(30:31)); [u1, u2] = tmp3{:};
    
    maxage = 101;
    totalcc1992 = y0(1);
    totaldc1992 = y0(2);
    totalhcc1992 = y0(3);

    x = (0:100)';
    v1=0.02*ones(101,1);
    v2=0.008*ones(101,1);
    v3=[0.009*ones(41,1); 0.021*ones(60,1)];
    v4=0.034*ones(101,1);
    v5=[0*ones(41,1); 0.003*ones(60,1)];
    m1x = @(x) m1 ./ (1 + exp(-m2 * (x - m3))) + m4;  dcmr = m1x(x);
    m2x = @(x) m5 ./ (1 + exp(-m6 * (x - m7))) + m8;  hccmr = m2x(x);
    p = [a0, a1, a1, a1, a1, a1, a2, a2, a2, a2, a2, a3, a3, a3, a3, a3, a4, a4, a4, a4, a4, a5, a5, a5, a5, a5, a6, a6, a6, a6, a6, a7, a7, a7, a7, a7, a8, a8, a8, a8, a8, a9, a9, a9, a9, a9, a10, a10, a10, a10, a10, a11, a11, a11, a11, a11, a12, a12, a12, a12, a12, a13, a13, a13, a13, a13, a14, a14, a14, a14, a14, a15, a15, a15, a15, a15, a16, a16, a16, a16, a16, a17, a17, a17, a17, a17, a18, a18, a18, a18, a18, a19, a19, a19, a19, a19, a20, a20, a20, a20, a20]';
    pdf_vals = normpdf(x, u1, u2); pdf_vals = pdf_vals(:);ageHCCshare = (pdf_vals) / sum(pdf_vals);

    chbpop = HBVdata20244(:,1) .* p;
    ccpop = HBVdata20241(1:101,15) * totalcc1992;
    dcpop = HBVdata20241(1:101,15) * totaldc1992;
    hccpop = ageHCCshare * totalhcc1992;

    zlthbvtmchb(:,1, m) = chbpop;  
    zlthbvtmcc(:,1, m) = ccpop;
    zlthbvtmdc(:,1, m) = dcpop; 
    zlthbvtmhcc(:,1, m) = hccpop; 

    for j = 1:N
        chbpop = HBVdata20244(:,j) .* p;
        totalpopa = HBVdata20245(:,j);
        totalpop = sum(totalpopa);
        tadr(:,j) = HBVdata20242(:,j);
        nonhbvmr(:,j) = (totalpopa .* 0.001 .* tadr(:,j) - dcpop .* dcmr - hccpop .* hccmr) ./ totalpopa;
        survival = 1 - nonhbvmr;
        ccnewpop(1) = 0; dcnewpop(1) = 0; hccnewpop(1) = 0;

        for i = 2:maxage
            ccnewpop(i)  = (survival(i-1,j) - v2(i-1) - v3(i-1)) * ccpop(i-1) + v1(i-1) * chbpop(i-1);
            dcnewpop(i)  = (survival(i-1,j) - v4(i-1) - dcmr(i-1)) * dcpop(i-1) + v2(i-1) * ccpop(i-1);
            hccnewpop(i) = (survival(i-1,j) - hccmr(i-1)) * hccpop(i-1) + v4(i-1) * dcpop(i-1) + v5(i-1) * chbpop(i-1) + v3(i-1) * ccpop(i-1);
            hccnew(i)    = v4(i-1) * dcpop(i-1) + v5(i-1) * chbpop(i-1) + v3(i-1) * ccpop(i-1);
            hccdeath(i)  = hccmr(i-1) * hccpop(i-1);
            dcdeath(i) = dcmr(i-1) * dcpop(i-1);
            dcnew(i)  = v2(i-1) * ccpop(i-1);
            ccnew(i) = v1(i-1) * chbpop(i-1);
        end

        for i = 1:maxage
            hccpop(i) = hccnewpop(i);
            ccpop(i)  = ccnewpop(i);
            dcpop(i)  = dcnewpop(i);
            hccdeathvec(i) = hccdeath(i);
            hccnewvec(i)   = hccnew(i);
            dcdeathvec(i) = dcdeath(i);
            dcnewvec(i)   = dcnew(i);
            ccnewvec(i) = ccnew(i);
        end

        zlthbvtmchb(:,j+1, m)   = chbpop;
        zlthbvtmcc(:,j+1, m)  = ccpop;
        zlthbvtmdc(:,j+1, m)  = dcpop; 
        zlthbvtmhcc(:,j+1, m) = hccpop; 
        hccdeathmat(:,j, m)   = hccdeathvec;
        hccnewmat(:,j, m)     = hccnewvec;
        dcdeathmat(:,j, m)   = dcdeathvec;
        dcnewmat(:,j, m)     = dcnewvec;
        ccnewmat(:,j, m)     = ccnewvec;
    end
    nonhbvmrmat(:,:,m) = nonhbvmr; 

    for j = 1:N-1
        hccmortality(1,j,m)  = sum(hccdeathmat(1,j,m))    / sum(HBVdata20245(1,j+1))    * 100000;
        hccmortality(2,j,m)  = sum(hccdeathmat(2:5,j,m))  / sum(HBVdata20245(2:5,j+1))  * 100000;
        hccmortality(3,j,m)  = sum(hccdeathmat(6:10,j,m)) / sum(HBVdata20245(6:10,j+1)) * 100000;
        hccmortality(4,j,m)  = sum(hccdeathmat(11:15,j,m)) / sum(HBVdata20245(11:15,j+1)) * 100000;
        hccmortality(5,j,m)  = sum(hccdeathmat(16:20,j,m)) / sum(HBVdata20245(16:20,j+1)) * 100000;
        hccmortality(6,j,m)  = sum(hccdeathmat(21:25,j,m)) / sum(HBVdata20245(21:25,j+1)) * 100000;
        hccmortality(7,j,m)  = sum(hccdeathmat(26:30,j,m)) / sum(HBVdata20245(26:30,j+1)) * 100000;
        hccmortality(8,j,m)  = sum(hccdeathmat(31:35,j,m)) / sum(HBVdata20245(31:35,j+1)) * 100000;
        hccmortality(9,j,m)  = sum(hccdeathmat(36:40,j,m)) / sum(HBVdata20245(36:40,j+1)) * 100000;
        hccmortality(10,j,m) = sum(hccdeathmat(41:45,j,m)) / sum(HBVdata20245(41:45,j+1)) * 100000;
        hccmortality(11,j,m) = sum(hccdeathmat(46:50,j,m)) / sum(HBVdata20245(46:50,j+1)) * 100000;
        hccmortality(12,j,m) = sum(hccdeathmat(51:55,j,m)) / sum(HBVdata20245(51:55,j+1)) * 100000;
        hccmortality(13,j,m) = sum(hccdeathmat(56:60,j,m)) / sum(HBVdata20245(56:60,j+1)) * 100000;
        hccmortality(14,j,m) = sum(hccdeathmat(61:65,j,m)) / sum(HBVdata20245(61:65,j+1)) * 100000;
        hccmortality(15,j,m) = sum(hccdeathmat(66:70,j,m)) / sum(HBVdata20245(66:70,j+1)) * 100000;
        hccmortality(16,j,m) = sum(hccdeathmat(71:75,j,m)) / sum(HBVdata20245(71:75,j+1)) * 100000;
        hccmortality(17,j,m) = sum(hccdeathmat(76:80,j,m)) / sum(HBVdata20245(76:80,j+1)) * 100000;
        hccmortality(18,j,m) = sum(hccdeathmat(81:85,j,m)) / sum(HBVdata20245(81:85,j+1)) * 100000;
        hccmortality(19,j,m) = sum(hccdeathmat(86:101,j,m)) / sum(HBVdata20245(86:101,j+1)) * 100000;
    end
end

for i = 1:19
    for j = 1:N-1
        samples = squeeze(hccmortality(i, j, :));  
        hccmortality01(i, j) = prctile(samples, 50);     
        hccmortality02(i, j) = prctile(samples, 2.5);  
        hccmortality03(i, j) = prctile(samples, 97.5); 
    end
end

years = 1:19;  
obs = HBVdata20241(1:19, 47:49);  
xticklabel = {'<1','1–4','5–9','10–14','15–19','20–24','25–29','30–34','35–39','40–44','45–49','50–54','55–59','60–64','65–69','70–74','75–79','80–84','≥85'};
figure('Position', [100, 100, 1400, 400]);
for j = 1:3
    subplot(1, 3, j); hold on;
    fill([years, fliplr(years)], [hccmortality02(:,25+j-1)', fliplr(hccmortality03(:,25+j-1)')], [0.8 0.8 1], 'EdgeColor', 'none');
    plot(years, hccmortality01(:,25+j-1), 'k-', 'LineWidth', 1.5);
    plot(years, obs(:,j), 'ro', 'MarkerFaceColor', 'r');
    title(num2str(2015 + j));
    xlabel('Age'); ylabel('HCC mortality (per 100,000)');
    xlim([1, 19]); ylim([0, 120]); yticks(0:20:120);
    xticks(1:19); xticklabels(xticklabel); xtickangle(45);
    set(gca, 'FontSize', 7); grid off;
end

names = {'zlthbvtmchb','zlthbvtmcc','zlthbvtmdc','zlthbvtmhcc'};
for idx = 1:length(names)
    name = names{idx};
    mat = eval(name);  
    [r, c, ~] = size(mat);  
    
    for i = 1:r
        for j = 1:c
            samples = squeeze(mat(i, j, :));
            m_val(i, j) = prctile(samples, 50);
            lower(i, j) = prctile(samples, 2.5);
            upper(i, j) = prctile(samples, 97.5);
        end
    end

    writematrix(m_val, fullfile(output_dir2, sprintf('1_%s101.csv', name)));
    writematrix(lower, fullfile(output_dir2, sprintf('1_%s102.csv', name)));
    writematrix(upper, fullfile(output_dir2, sprintf('1_%s103.csv', name)));
end

names = ["hccdeathmat", "hccnewmat","dcdeathmat", "dcnewmat","nonhbvmrmat","ccnewmat"];
for idx = 1:length(names)
    name = names{idx};
    mat = eval(name);  
    [r, c, ~] = size(mat);  
    
    for i = 1:r
        for j = 1:c
            samples = squeeze(mat(i, j, :));
            m_val(i, j) = prctile(samples, 50);
            lower(i, j) = prctile(samples, 2.5);
            upper(i, j) = prctile(samples, 97.5);
        end
    end

    writematrix(m_val, fullfile(output_dir2, sprintf('1_%s201.csv', name)));
    writematrix(lower, fullfile(output_dir2, sprintf('1_%s202.csv', name)));
    writematrix(upper, fullfile(output_dir2, sprintf('1_%s203.csv', name)));
end

tmp1 = num2cell(theta_vals(22:29)); [m1, m2, m3, m4, m5, m6, m7, m8] = tmp1{:};
m1x = @(x) m1 ./ (1 + exp(-m2 * (x - m3))) + m4;  dcmr = m1x(x);
m2x = @(x) m5 ./ (1 + exp(-m6 * (x - m7))) + m8;  hccmr = m2x(x);

save_filename = fullfile(output_dir2, 'HBV1.mat');
save(save_filename, 'zlthbvtmchb', 'zlthbvtmcc', 'zlthbvtmdc', 'zlthbvtmhcc', ...
                    'hccdeathmat', 'hccnewmat', 'dcdeathmat', 'dcnewmat', ...
                    'nonhbvmrmat', 'ccnewmat', '-v6');