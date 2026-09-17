function hccmortality = mcmc_model(theta, y0, data, N)

global zlthbvtmchb zlthbvtmdc zlthbvtmhcc hccdeathmat hccnewmat zlthbvtmcc

HBVdata20241 = data.HBVdata20241;
HBVdata20242 = data.HBVdata20242;
HBVdata20243 = data.HBVdata20243;
HBVdata20244 = data.HBVdata20244;
HBVdata20245 = data.HBVdata20245;

tmp1 = num2cell(theta(22:29)); [m1, m2, m3, m4, m5, m6, m7, m8] = tmp1{:};
tmp2 = num2cell(theta(1:21)); [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20] = tmp2{:};
tmp3 = num2cell(theta(30:31)); [u1, u2] = tmp3{:};

maxage=101;
totalcc1992=y0(1);
totaldc1992=y0(2);
totalhcc1992=y0(3);

x=(0:100)';
v1=0.02*ones(101,1);
v2=0.008*ones(101,1);
v3=[0.009*ones(41,1); 0.021*ones(60,1)];
v4=0.034*ones(101,1);
v5=[0*ones(41,1); 0.003*ones(60,1)];
m1x = @(x) m1 ./ (1 + exp(-m2 * (x - m3))) + m4;  dcmr = m1x(x);
m2x = @(x) m5 ./ (1 + exp(-m6 * (x - m7))) + m8;  hccmr = m2x(x);
p = [a0, a1, a1, a1, a1, a1, a2, a2, a2, a2, a2, a3, a3, a3, a3, a3, a4, a4, a4, a4, a4, a5, a5, a5, a5, a5, a6, a6, a6, a6, a6, a7, a7, a7, a7, a7, a8, a8, a8, a8, a8, a9, a9, a9, a9, a9, a10, a10, a10, a10, a10, a11, a11, a11, a11, a11, a12, a12, a12, a12, a12, a13, a13, a13, a13, a13, a14, a14, a14, a14, a14, a15, a15, a15, a15, a15, a16, a16, a16, a16, a16, a17, a17, a17, a17, a17, a18, a18, a18, a18, a18, a19, a19, a19, a19, a19, a20, a20, a20, a20, a20]';
pdf_vals = normpdf(x, u1, u2); pdf_vals = pdf_vals(:);ageHCCshare = (pdf_vals) / sum(pdf_vals);

chbpop=HBVdata20244(:,1).* p;
ccpop=HBVdata20241(1:101,15)*totalcc1992;
dcpop=HBVdata20241(1:101,15)*totaldc1992;
hccpop = ageHCCshare*totalhcc1992;

zlthbvtmchb(:,1)=chbpop;  
zlthbvtmcc(:,1)=ccpop;
zlthbvtmdc(:,1)=dcpop; 
zlthbvtmhcc(:,1)=hccpop; 

for j=1:N
    chbpop=HBVdata20244(:,j).* p;
    totalpopa=HBVdata20245(:,j);
    totalpop=sum(totalpopa);
    tadr(:,j)=HBVdata20242(:,j); 
    nonhbvmr(:,j)=(totalpopa.*0.001.*tadr(:,j)-dcpop.*dcmr-hccpop.*hccmr)./totalpopa;

    survival=1-nonhbvmr;
    ccnewpop(1)=0;
    dcnewpop(1)=0;
    hccnewpop(1)=0;
    for i=2:maxage
        ccnewpop(i)=(survival(i-1,j)-v2(i-1)-v3(i-1))*ccpop(i-1)+v1(i-1)*chbpop(i-1);
        dcnewpop(i)=(survival(i-1,j)-v4(i-1)-dcmr(i-1))*dcpop(i-1)+v2(i-1)*ccpop(i-1);
        hccnewpop(i)=(survival(i-1,j)-hccmr(i-1))*hccpop(i-1)+v4(i-1)*dcpop(i-1)+v5(i-1)*chbpop(i-1)+v3(i-1)*ccpop(i-1);
        hccdeath(i)=(1-survival(i-1,j)+hccmr(i-1))*hccpop(i-1);
        hccnew(i)=v4(i-1)*dcpop(i-1)+v5(i-1)*chbpop(i-1)+v3(i-1)*ccpop(i-1);
    end
    for i=1:maxage
        hccpop(i)=hccnewpop(i);
        ccpop(i)=ccnewpop(i);
        dcpop(i)=dcnewpop(i);
        hccdeathvec(i)=hccdeath(i);
        hccnewvec(i)=hccnew(i);
    end
    zlthbvtmchb(:,j)=chbpop;  
    zlthbvtmcc(:,j+1)=ccpop;
    zlthbvtmdc(:,j+1)=dcpop; 
    zlthbvtmhcc(:,j+1)=hccpop; 
    hccdeathmat(:,j)=hccdeathvec;
    hccnewmat(:,j)=hccnewvec;
end
for j=1:N
    hccmortality(1,j)  = sum(hccdeathmat(1,j))   / sum(HBVdata20245(1,j+1))    * 100000;
    hccmortality(2,j)  = sum(hccdeathmat(2:5,j)) / sum(HBVdata20245(2:5,j+1))  * 100000;
    hccmortality(3,j)  = sum(hccdeathmat(6:10,j)) / sum(HBVdata20245(6:10,j+1)) * 100000;
    hccmortality(4,j)  = sum(hccdeathmat(11:15,j)) / sum(HBVdata20245(11:15,j+1)) * 100000;
    hccmortality(5,j)  = sum(hccdeathmat(16:20,j)) / sum(HBVdata20245(16:20,j+1)) * 100000;
    hccmortality(6,j)  = sum(hccdeathmat(21:25,j)) / sum(HBVdata20245(21:25,j+1)) * 100000;
    hccmortality(7,j)  = sum(hccdeathmat(26:30,j)) / sum(HBVdata20245(26:30,j+1)) * 100000;
    hccmortality(8,j)  = sum(hccdeathmat(31:35,j)) / sum(HBVdata20245(31:35,j+1)) * 100000;
    hccmortality(9,j)  = sum(hccdeathmat(36:40,j)) / sum(HBVdata20245(36:40,j+1)) * 100000;
    hccmortality(10,j) = sum(hccdeathmat(41:45,j)) / sum(HBVdata20245(41:45,j+1)) * 100000;
    hccmortality(11,j) = sum(hccdeathmat(46:50,j)) / sum(HBVdata20245(46:50,j+1)) * 100000;
    hccmortality(12,j) = sum(hccdeathmat(51:55,j)) / sum(HBVdata20245(51:55,j+1)) * 100000;
    hccmortality(13,j) = sum(hccdeathmat(56:60,j)) / sum(HBVdata20245(56:60,j+1)) * 100000;
    hccmortality(14,j) = sum(hccdeathmat(61:65,j)) / sum(HBVdata20245(61:65,j+1)) * 100000;
    hccmortality(15,j) = sum(hccdeathmat(66:70,j)) / sum(HBVdata20245(66:70,j+1)) * 100000;
    hccmortality(16,j) = sum(hccdeathmat(71:75,j)) / sum(HBVdata20245(71:75,j+1)) * 100000;
    hccmortality(17,j) = sum(hccdeathmat(76:80,j)) / sum(HBVdata20245(76:80,j+1)) * 100000;
    hccmortality(18,j) = sum(hccdeathmat(81:85,j)) / sum(HBVdata20245(81:85,j+1)) * 100000;
    hccmortality(19,j) = sum(hccdeathmat(86:101,j)) / sum(HBVdata20245(86:101,j+1)) * 100000;
end
hccmortality = hccmortality(:,13:23);
end