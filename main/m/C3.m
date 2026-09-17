clear all;
clear model data params options
rng(888);

input_dir  = fullfile('..', 'datasets');
output_dir = fullfile('..', 'results', 'C1');

N = 133;
maxage = 101;
nMC = 5000;   

rc_ci = [
    0.0801  0.0144  0.1458
    0.0462  0.0120  0.0804
    0.0056  0.0000  0.0373
    0.0144  0.0000  0.0487
    0.0437  0.0152  0.0721
    0.0000  0.0000  0.0224
    0.0000  0.0000  0.0266
    0.0342  0.0090  0.0595
    0.0000  0.0000  0.0107
    0.0258  0.0113  0.0403
    0.0180  0.0154  0.0206
];

rc_lower = rc_ci(:,2);
rc_upper = rc_ci(:,3);
num_rc = size(rc_ci,1);
rc_samples = zeros(num_rc, nMC);
recoverrate_samples = zeros(maxage, nMC);

for i = 1:nMC
    rc = rc_lower + (rc_upper - rc_lower) .* rand(num_rc,1);
    rc1 = rc(1);  rc2 = rc(2);  rc3 = rc(3);  rc4 = rc(4);  rc5 = rc(5);
    rc6 = rc(6);  rc7 = rc(7);  rc8 = rc(8);  rc9 = rc(9);  rc10 = rc(10);
    rc11 = rc(11);

    recoverrate_samples(:,i) = [rc11,repmat(rc1,1,4), repmat(rc2,1,5),repmat(rc3,1,5), repmat(rc4,1,5),...
        repmat(rc5,1,5), repmat(rc6,1,5),repmat(rc7,1,5),  repmat(rc8,1,5),repmat(rc9,1,10), repmat(rc10,1,10),repmat(rc11,1,41)]';
end

m2 = matfile(fullfile(input_dir,'data1.mat'));
data = m2.data1;
transmrate_samples = zeros(maxage, nMC);
HBVdata20241 = data.HBVdata20241;
lower_bound = HBVdata20241(1:101,13); 
upper_bound = HBVdata20241(1:101,14);  
for mc = 1:nMC
    transmrate_samples(:, mc) = lower_bound + (upper_bound - lower_bound) .* rand(maxage, 1);
end

m1 = matfile(fullfile(input_dir,'data1.mat'));
data = m1.data1;
HBVdata20241 = data.HBVdata20241;
HBVdata20242 = data.HBVdata20242;
HBVdata20243 = data.HBVdata20243;
birthrate = HBVdata20241(1:N,1); 
intranterine = 0.06;  
vaccinarate = 0.85 * HBVdata20241(1:N,2) ; 
disedeath     = HBVdata20241(1:maxage,11);
sctransmrate  = HBVdata20241(1:maxage,4);

zlthbvtm1 = zeros(maxage, N+1, nMC);
zlthbvtm2 = zeros(maxage, N+1, nMC);
zlthbvtm3 = zeros(maxage, N+1, nMC);
zlthbvtm4 = zeros(maxage, N,   nMC);
zlthbvtm5 = zeros(maxage, N,   nMC);
totppozu = zeros(maxage, N+1, nMC);
HBVdeath_all = zeros(maxage, N, nMC);
yrhbsrate_all = zeros(N, 1, nMC);

for m = 1:nMC
    
    recoverrate1 = recoverrate_samples(:, m);
    for year=1:N
        recoverrate(:,year)=1.0.*recoverrate1; 
    end

    transmrate(:,1)=transmrate_samples(:,m);
    rc131=0.5308; rc14=0.3853;
    for y=1:N
        transmrate(:,y+1)=(rc131*exp(-rc14*y)+(1-rc131)).*transmrate(:,1);
    end

    spop = HBVdata20241(1:maxage,6);
    cpop = HBVdata20241(1:maxage,7);
    rpop = HBVdata20241(1:maxage,8);

    zlthbvtm1(:,1, m) = spop;
    zlthbvtm2(:,1, m) = cpop;
    zlthbvtm3(:,1, m) = rpop;

    totppozu(:,1, m)=zlthbvtm1(:,1, m)+zlthbvtm2(:,1, m)+zlthbvtm3(:,1, m);

    for j=1:N
        totalpop=sum(spop)+sum(cpop)+sum(rpop);
        totalpopa=spop+cpop+rpop;
        totalcpop=sum(cpop);
        nadeath(:,j)=HBVdata20242(:,j);
        adultvaccinarate(:,j)=0.85.*HBVdata20243(:,j);  
        dedisedeath(:,j)=(totalpopa.*0.001.*nadeath(:,j)-cpop.*disedeath)./totalpopa;
        survival=1-dedisedeath;
        womanhbvrate(j)=intranterine*0.8491*sum(cpop(16:50))/sum(totalpopa(16:50));
        snewpop(1)=birthrate(j)*totalpop*(1-vaccinarate(j)-womanhbvrate(j));
        cnewpop(1)=birthrate(j)*totalpop*0.9*womanhbvrate(j);
        rnewpop(1)=birthrate(j)*totalpop*vaccinarate(j)+birthrate(j)*totalpop*0.1*womanhbvrate(j);
        newnum(j)=birthrate(j)*totalpop;
        HBVdeath(1)=0;
        newinfnum(1)=birthrate(j)*totalpop*0.9*womanhbvrate(j);

        for i=2:maxage
            snewpop(i)=survival(i-1,j)*spop(i-1) -adultvaccinarate(i-1,j)*spop(i-1) -transmrate(i-1,j+1)*totalcpop*spop(i-1)/totalpop;
            cnewpop(i)=(survival(i-1,j)-recoverrate(i-1,j)-disedeath(i-1))*cpop(i-1) +sctransmrate(i-1)*transmrate(i-1,j+1)*totalcpop*spop(i-1)/totalpop;
            rnewpop(i)=survival(i-1,j)*rpop(i-1) +recoverrate(i-1,j)*cpop(i-1) +adultvaccinarate(i-1,j)*spop(i-1) +(1-sctransmrate(i-1))*transmrate(i-1,j+1)*totalcpop*spop(i-1)/totalpop;
            HBVdeath(i)=disedeath(i-1)*cpop(i-1);
            newinfnum(i)=sctransmrate(i-1)*transmrate(i-1,j+1)*totalcpop*(spop(i-1))/(totalpop);  
        end

        for i=1:maxage
            spop(i)=snewpop(i);
            cpop(i)=cnewpop(i);
            rpop(i)=rnewpop(i);
            hbsagrate(i)=100*cpop(i)/(spop(i)+cpop(i)+rpop(i));
            HBVdeatha(i)=HBVdeath(i);
            newc(i)=newinfnum(i);
        end

       zlthbvtm1(:,j+1, m)   = spop; 
       zlthbvtm2(:,j+1, m)   = cpop;
       zlthbvtm3(:,j+1, m)   = rpop;
       zlthbvtm4(:,j, m)  = totalpopa;
       HBVdeath_all(:,j, m)  = HBVdeatha;
       zlthbvtm5(:,j, m)  = newc;
       totppozu(:,j+1, m)=zlthbvtm2(:,j+1, m)+zlthbvtm1(:,j+1, m)+zlthbvtm3(:,j+1, m);
       yrhbsrate_all(1,1,m) = 100 * sum(zlthbvtm2(2:60,1, m)) / sum(totppozu(2:60,1, m));
       for j = 1:N-1
           yrhbsrate_all(j+1,1,m) = 100 * sum(zlthbvtm2(2:60,j+1, m)) / sum(totppozu(2:60,j+1, m));
       end
    end
end

a = [0.001323959, 0.327029632, 0.045000807, 0.044706294, 0.405412702, ...
     0.684596414, 0.992912748, 0.079065542, 0.885114999, 0.492711747, ...
     0.943946583, 0.039305817, 0.359766887, 0.00506399, 0.002765995, ...
     0.682467642, 0.74539068, 0.202819098, 0.357055858, 0.114390501, ...
     0.165575045];
p = [a(1), repelem(a(2:21), 5)]'; 

zlthbvtm6 = zeros(maxage, N, nMC);
for m = 1:nMC
    zlthbvtm6(:, :, m) = zlthbvtm5(:, :, m) .* p;
end

names = {'zlthbvtm2','zlthbvtm4','HBVdeath_all',"zlthbvtm5"};
for idx = 1:length(names)
    name = names{idx};
    mat = eval(name);  
    [r, c, ~] = size(mat);  
    
    for i = 1:r
        for j = 1:c
            samples = squeeze(mat(i, j, :));
            m(i, j) =  prctile(samples, 50);
            lower(i, j) = prctile(samples, 2.5);
            upper(i, j) = prctile(samples, 97.5);
        end
    end

    writematrix(m,     fullfile(output_dir, sprintf('3_%s101.csv', name)));
    writematrix(lower, fullfile(output_dir, sprintf('3_%s102.csv', name)));
    writematrix(upper, fullfile(output_dir, sprintf('3_%s103.csv', name)));
end

mat = yrhbsrate_all;  
[r, c, nMC] = size(mat);  
m = zeros(r,c);
lower = zeros(r,c);
upper = zeros(r,c);
for i = 1:r
    for j = 1:c  
        samples = squeeze(mat(i, j, :));  
        m(i,j) =  prctile(samples, 50);
        lower(i,j) = prctile(samples, 2.5);
        upper(i,j) = prctile(samples, 97.5);
    end
end
writematrix(m,     fullfile(output_dir, sprintf('3_%s101.csv', "yrhbsrate")));
writematrix(lower, fullfile(output_dir, sprintf('3_%s102.csv', "yrhbsrate")));
writematrix(upper, fullfile(output_dir, sprintf('3_%s103.csv', "yrhbsrate")));
save_filename = fullfile(output_dir, 'C3.mat');
save(save_filename, 'zlthbvtm2', 'zlthbvtm4', 'HBVdeath_all', 'zlthbvtm5','zlthbvtm6', 'yrhbsrate_all', '-v6');