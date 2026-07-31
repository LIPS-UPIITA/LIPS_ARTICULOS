% ConcordanceMetrics.m
% Reproduce Cohen's kappa and Gwet's AC1 for the provided 2x2 table
% Verification with Bootstrap method for small sample sizes.

% Author: LIGJ,26/2026
% v1.2
clear; clc;

% Contingency counts
a = 11; % R1=1, R2=1
b = 1;  % R1=1, R2=0
c = 5;  % R1=0, R2=1
d = 0;  % R1=0, R2=0
N = a + b + c + d;

% Observed agreement
Po = (a + d) / N;

% Marginal proportions
p1_yes = (a + b) / N;
p1_no  = (c + d) / N;
p2_yes = (a + c) / N;
p2_no  = (b + d) / N;

% Cohen's kappa
Pe_kappa = p1_yes * p2_yes + p1_no * p2_no;
kappa = (Po - Pe_kappa) / (1 - Pe_kappa);

% Approximate SE for kappa (the simplified approximation used in the doc)
SE_kappa = sqrt( Po * (1 - Po) / ( N * (1 - Pe_kappa)^2 ) );

z_kappa = kappa / SE_kappa;
p_kappa = 2 * (1 - normcdf(abs(z_kappa))); % two-sided
CI95_kappa = [kappa - 1.96*SE_kappa, kappa + 1.96*SE_kappa];

% Gwet's AC1
pbar_yes = (p1_yes + p2_yes) / 2;
pbar_no  = (p1_no  + p2_no)  / 2;
Pe_AC1 = pbar_yes^2 + pbar_no^2;
AC1 = (Po - Pe_AC1) / (1 - Pe_AC1);

% Display results
fprintf('N = %d\n', N);
fprintf('Po = %.8f\n\n', Po);

fprintf('Cohen''s kappa:\n');
fprintf(' Pe_kappa = %.8f\n', Pe_kappa);
fprintf(' kappa = %.8f\n', kappa);
fprintf(' SE_kappa (approx) = %.8f\n', SE_kappa);
fprintf(' z = %.4f, p = %.4f\n', z_kappa, p_kappa);
fprintf(' 95%% CI (approx) = [%.4f, %.4f]\n\n', CI95_kappa(1), CI95_kappa(2));

fprintf('Gwet''s AC1:\n');
fprintf(' pbar_yes = %.8f, pbar_no = %.8f\n', pbar_yes, pbar_no);
fprintf(' Pe_AC1 = %.8f\n', Pe_AC1);
fprintf(' AC1 = %.8f\n\n', AC1);

% Verificaction: bootstrap CI for kappa and AC1, mainly recommended for small N
nBoot = 20000;
idxAll = repmat([1:N], nBoot, 1);
rng(0);
bootK = zeros(nBoot,1);
bootAC1 = zeros(nBoot,1);
counts = [a,b,c,d];
% Create vector of pairs to resample: [R1,R2] for each item
pairs = [ones(a,1)*1, ones(a,1)*1;
         ones(b,1)*1, zeros(b,1)*1;
         zeros(c,1)*1, ones(c,1)*1;
         zeros(d,1)*1, zeros(d,1)*1];
for i=1:nBoot
    s = randsample(N,N,true);
    samp = pairs(s,:);
    a_s = sum(samp(:,1)==1 & samp(:,2)==1);
    b_s = sum(samp(:,1)==1 & samp(:,2)==0);
    c_s = sum(samp(:,1)==0 & samp(:,2)==1);
    d_s = sum(samp(:,1)==0 & samp(:,2)==0);
    Po_s = (a_s + d_s)/N;
    p1_yes_s = (a_s + b_s)/N; p1_no_s = (c_s + d_s)/N;
    p2_yes_s = (a_s + c_s)/N; p2_no_s = (b_s + d_s)/N;
    Pe_k_s = p1_yes_s * p2_yes_s + p1_no_s * p2_no_s;
    bootK(i) = (Po_s - Pe_k_s) / (1 - Pe_k_s);
    pbar_yes_s = (p1_yes_s + p2_yes_s)/2;
    pbar_no_s = (p1_no_s + p2_no_s)/2;
    Pe_ac1_s = pbar_yes_s^2 + pbar_no_s^2;
    bootAC1(i) = (Po_s - Pe_ac1_s) / (1 - Pe_ac1_s);
end
CI_kappa_boot = prctile(bootK, [2.5 97.5]);
CI_AC1_boot = prctile(bootAC1, [2.5 97.5]);
fprintf('Bootstrap CI kappa (95%%): [%.4f, %.4f]\n', CI_kappa_boot);
fprintf('Bootstrap CI AC1 (95%%): [%.4f, %.4f]\n', CI_AC1_boot);
