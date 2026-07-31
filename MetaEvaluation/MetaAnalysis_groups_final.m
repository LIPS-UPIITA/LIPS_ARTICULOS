%% MetaAnalysis_groups_final.m
% MATLAB R2024b compatible
%
% Meta-analysis of Mean Differences (MD)
% Random-effects models:
%   1) DerSimonian-Laird (DL)
%   2) REML
%
% Three prespecified subgroup-specific meta-analyses:
%   - Elderly_Female: k = 2
%   - Elderly_Male:   k = 2
%   - Elderly:        k = 4
%
% IMPORTANT:
% No overall pooled estimate is calculated across the 8 rows because
% some studies contribute to more than one subgroup analysis.
%
% Prediction intervals:
%   - Calculated only when k >= 3
%   - Not calculated for k < 3
%   - Results with k < 5 should be interpreted cautiously.
%
% Sensitivity analysis:
%   - Exclude Hu(2021) from Elderly group.
%
% Effect direction:
%   MD = Mean_Control - Mean_Sarcopenia
%   Positive MD = higher value in controls than in sarcopenia.
%
% MATLAB R2024b
% -------------------------------------------------------------------------

clear;
close all;
clc;

%% ========================================================================
% DATA
% Each row:
% {StudyName, GroupLabel, nControl, meanControl, sdControl,
%  nCase, meanCase, sdCase}
% ========================================================================

T = {
'He(2024)_F',     'Elderly_Female', 68, 20.4, 3.6, 36, 14.8, 2.6;
'Li(2024)_F',     'Elderly_Female', 29, 22.4, 3.0, 31, 14.4, 2.8;

'He(2024)_M',     'Elderly_Male',   30, 32.1, 5.3, 19, 24.1, 2.9;
'Li(2024)_M',     'Elderly_Male',   16, 34.1, 5.4, 17, 23.8, 2.2;

'Sepulveda(2025)','Elderly',         22, 24.5, 7.4, 13, 19.8, 5.9;
'He(2024)',       'Elderly',         98, 23.9, 7.0, 55, 18.0, 5.2;
'Li(2024)',       'Elderly',         45, 26.6, 6.9, 48, 17.7, 5.2;
'Hu(2021)',       'Elderly',          5, 29.6, 8.88, 5, 20.0, 3.74;
};

%% ========================================================================
% CONVERT DATA
% ========================================================================

studies = T(:,1);
groups  = T(:,2);

nC  = cell2mat(T(:,3));
mC  = cell2mat(T(:,4));
sdC = cell2mat(T(:,5));

nT  = cell2mat(T(:,6));
mT  = cell2mat(T(:,7));
sdT = cell2mat(T(:,8));

%% ========================================================================
% STUDY-LEVEL EFFECT SIZES
% ========================================================================

% Mean Difference:
% Control - Sarcopenia
MD = mC - mT;

% Within-study variance of the MD
var_within = (sdC.^2)./nC + (sdT.^2)./nT;

% Standard error
se_within = sqrt(var_within);

%% ========================================================================
% OUTPUT DIRECTORY
% ========================================================================

outdir = fullfile(pwd,'outputs');

if ~exist(outdir,'dir')
    mkdir(outdir);
end

%% ========================================================================
% STUDY-LEVEL TABLE
% ========================================================================

T_out = table( ...
    studies(:), ...
    groups(:), ...
    nC(:), ...
    nT(:), ...
    MD(:), ...
    se_within(:), ...
    var_within(:), ...
    'VariableNames', ...
    {'Study','Subgroup','nControl','nCase','MD','SE','Var'});

writetable( ...
    T_out, ...
    fullfile(outdir,'per_study_table_allrows.xlsx'));

%% ========================================================================
% PROCESS EACH PREDEFINED GROUP
% ========================================================================

uniqueGroups = unique(groups,'stable');

OutGroups = struct();

summary_rows = {};

for ig = 1:numel(uniqueGroups)

    grp = uniqueGroups{ig};

    fprintf('\n============================================================\n');
    fprintf('Processing group: %s\n',grp);
    fprintf('============================================================\n');

    %% ---------------------------------------------------------------
    % Select studies belonging to this subgroup
    % ---------------------------------------------------------------

    idx = find(strcmp(groups,grp));

    y = MD(idx)';
    v = var_within(idx)';

    k = numel(y);

    fprintf('Number of studies: k = %d\n',k);

    if k == 0
        warning('No studies found for group %s. Skipping.',grp);
        continue;
    end

    %% ---------------------------------------------------------------
    % FIXED-EFFECT WEIGHTS
    % ---------------------------------------------------------------

    w_fixed = 1 ./ v;

    theta_fixed = ...
        sum(w_fixed .* y) / sum(w_fixed);

    %% ---------------------------------------------------------------
    % COCHRAN Q
    % ---------------------------------------------------------------

    Q = sum( ...
        w_fixed .* ...
        (y - theta_fixed).^2);

    df_Q = k - 1;

    %% ---------------------------------------------------------------
    % DL CONSTANT C
    % ---------------------------------------------------------------

    C = ...
        sum(w_fixed) - ...
        (sum(w_fixed.^2) / sum(w_fixed));

    %% ---------------------------------------------------------------
    % DER SIMONIAN-LAIRD TAU^2
    % ---------------------------------------------------------------

    if k > 1 && C > 0

        tau2_DL = max( ...
            0, ...
            (Q - df_Q) / C);

    else

        tau2_DL = 0;

    end

    %% ---------------------------------------------------------------
    % DL POOLED EFFECT
    % ---------------------------------------------------------------

    [theta_DL,var_theta_DL] = ...
        pooled_given_tau(y,v,tau2_DL);

    se_theta_DL = sqrt(var_theta_DL);

    ci_DL = ...
        theta_DL + ...
        norminv([0.025 0.975]) .* se_theta_DL;

    %% ---------------------------------------------------------------
    % DL WEIGHTS
    % ---------------------------------------------------------------

    wi_DL = 1 ./ (v + tau2_DL);

    weights_pct_DL = ...
        100 .* wi_DL ./ sum(wi_DL);

    %% ---------------------------------------------------------------
    % I2
    % ---------------------------------------------------------------

    if Q > 0 && k > 1

        I2 = max( ...
            0, ...
            ((Q - df_Q) / Q) * 100);

    else

        I2 = 0;

    end

    %% =================================================================
    % REML ESTIMATION
    % =================================================================
    % ============================================================
% REML estimation of tau^2
% ============================================================

opt = optimset( ...
    'TolX', 1e-10, ...
    'MaxIter', 1000, ...
    'Display', 'off');

% Conservative upper bound for tau^2
upper = max([ ...
    100 * max(v), ...
    100 * var(y), ...
    1]);

fun = @(t) restrictedLogLik(t, y, v);

try

    tau2_REML = fminbnd(fun, 0, upper, opt);

    % Numerical protection
    tau2_REML = max(0, tau2_REML);

catch ME

    warning( ...
        'REML optimization failed for group %s: %s. tau2_REML set to 0.', ...
        grp, ME.message);

    tau2_REML = 0;

end
    
    %% ---------------------------------------------------------------
    % REML POOLED EFFECT
    % ---------------------------------------------------------------

    [theta_REML,var_theta_REML] = ...
        pooled_given_tau( ...
        y,v,tau2_REML);

    se_theta_REML = sqrt(var_theta_REML);

    ci_REML = ...
        theta_REML + ...
        norminv([0.025 0.975]) .* se_theta_REML;

    %% ---------------------------------------------------------------
    % REML WEIGHTS
    % ---------------------------------------------------------------

    wi_REML = ...
        1 ./ (v + tau2_REML);

    weights_pct_REML = ...
        100 .* wi_REML ./ sum(wi_REML);

    %% =================================================================
    % PREDICTION INTERVALS
    % =================================================================
    %
    % PI is NOT calculated for k < 3.
    %
    % For k >= 3:
    %
    % PI = theta +/- t_(k-2,0.975)
    %      * sqrt(SE_theta^2 + tau^2)
    %
    % NOTE:
    % Prediction intervals with very small k remain highly uncertain.
    % A warning is generated when k < 5.
    % =================================================================

    PI_DL = [NaN NaN];
    PI_REML = [NaN NaN];

    PI_status = "Not calculated";

    if k >= 3

        if k < 5

            warning( ...
                ['Prediction interval for group %s is based on only ' ...
                 'k=%d studies and should be interpreted cautiously.'], ...
                 grp,k);

            PI_status = "Calculated with k<5; interpret cautiously";

        else

            PI_status = "Calculated";

        end

        df_PI = k - 2;

        tcrit = tinv(0.975,df_PI);

        %% DL PI

        PI_se_DL = ...
            sqrt(var_theta_DL + tau2_DL);

        PI_DL = ...
            theta_DL + ...
            [-1 1] .* tcrit .* PI_se_DL;

        %% REML PI

        PI_se_REML = ...
            sqrt(var_theta_REML + tau2_REML);

        PI_REML = ...
            theta_REML + ...
            [-1 1] .* tcrit .* PI_se_REML;

    end

    %% =================================================================
    % STORE GROUP RESULTS
    % =================================================================

    G = struct();

    G.group = grp;
    G.k = k;
    G.idx = idx;

    G.studies = studies(idx);

    G.y = y;
    G.v = v;

    G.Q = Q;
    G.df_Q = df_Q;
    G.I2 = I2;

    % DL
    G.tau2_DL = tau2_DL;
    G.theta_DL = theta_DL;
    G.se_theta_DL = se_theta_DL;
    G.ci_DL = ci_DL;
    G.weights_pct_DL = weights_pct_DL;

    % REML
    G.tau2_REML = tau2_REML;
    G.theta_REML = theta_REML;
    G.se_theta_REML = se_theta_REML;
    G.ci_REML = ci_REML;
    G.weights_pct_REML = weights_pct_REML;

    % Prediction intervals
    G.PI_DL = PI_DL;
    G.PI_REML = PI_REML;
    G.PI_status = PI_status;

    %% Store structure

    fieldName = ...
        matlab.lang.makeValidName(grp);

    OutGroups.(fieldName) = G;

    %% =================================================================
    % SUMMARY TABLE
    % =================================================================

    summary_rows(end+1,:) = { ...
        grp, ...
        k, ...
        theta_DL, ...
        ci_DL(1), ...
        ci_DL(2), ...
        PI_DL(1), ...
        PI_DL(2), ...
        tau2_DL, ...
        theta_REML, ...
        ci_REML(1), ...
        ci_REML(2), ...
        PI_REML(1), ...
        PI_REML(2), ...
        tau2_REML, ...
        Q, ...
        df_Q, ...
        I2, ...
        PI_status};

    %% =================================================================
    % FOREST PLOTS
    % =================================================================

    makeForestGroup( ...
        'DL', ...
        G, ...
        outdir);

    makeForestGroup( ...
        'REML', ...
        G, ...
        outdir);

    %% =================================================================
    % CONSOLE OUTPUT
    % =================================================================

    fprintf('\nGroup: %s\n',grp);
    fprintf('k = %d\n',k);

    fprintf( ...
        'DL   : theta = %.3f [%.3f, %.3f], tau2 = %.4f\n', ...
        theta_DL, ...
        ci_DL(1), ...
        ci_DL(2), ...
        tau2_DL);

    fprintf( ...
        'REML : theta = %.3f [%.3f, %.3f], tau2 = %.4f\n', ...
        theta_REML, ...
        ci_REML(1), ...
        ci_REML(2), ...
        tau2_REML);

    fprintf( ...
        'Q = %.3f, df = %d, I2 = %.2f%%\n', ...
        Q,df_Q,I2);

    if k >= 3

        fprintf( ...
            'DL PI   = [%.3f, %.3f]\n', ...
            PI_DL(1),PI_DL(2));

        fprintf( ...
            'REML PI = [%.3f, %.3f]\n', ...
            PI_REML(1),PI_REML(2));

    else

        fprintf( ...
            'Prediction intervals: NOT calculated (k < 3)\n');

    end

end

%% ========================================================================
% SUMMARY TABLE
% ========================================================================

Tsum = cell2table( ...
    summary_rows, ...
    'VariableNames',{ ...
    'Group', ...
    'k', ...
    'Theta_DL', ...
    'CIlo_DL', ...
    'CIhi_DL', ...
    'PIlo_DL', ...
    'PIhi_DL', ...
    'tau2_DL', ...
    'Theta_REML', ...
    'CIlo_REML', ...
    'CIhi_REML', ...
    'PIlo_REML', ...
    'PIhi_REML', ...
    'tau2_REML', ...
    'Q', ...
    'df', ...
    'I2', ...
    'PI_status'});

writetable( ...
    Tsum, ...
    fullfile(outdir,'summary_by_group.xlsx'));

%% ========================================================================
% SAVE MATLAB RESULTS
% ========================================================================

save( ...
    fullfile(outdir, ...
    'MetaAnalysis_groups_results.mat'), ...
    'OutGroups', ...
    'T_out', ...
    'Tsum');

%% ========================================================================
% SENSITIVITY ANALYSIS
% EXCLUDE HU(2021) FROM ELDERLY GROUP
% ========================================================================

groupName = 'Elderly';

fieldName = ...
    matlab.lang.makeValidName(groupName);

if ~isfield(OutGroups,fieldName)

    warning( ...
        'Group %s not found. Sensitivity analysis skipped.', ...
        groupName);

else

    S = OutGroups.(fieldName);

    idx = S.idx;

    %% ---------------------------------------------------------------
    % Locate Hu(2021)
    % ---------------------------------------------------------------

    hu_label = 'Hu(2021)';

    hu_global_idx = ...
        find(strcmp(studies,hu_label),1);

    if isempty(hu_global_idx)

        warning( ...
            'Hu(2021) was not found.');

    elseif ~ismember(hu_global_idx,idx)

        warning( ...
            'Hu(2021) is not part of the Elderly group.');

    else

        %% -----------------------------------------------------------
        % Dataset without Hu
        % -----------------------------------------------------------

        idx_no_hu = ...
            idx(idx ~= hu_global_idx);

        y_s = MD(idx_no_hu)';
        v_s = var_within(idx_no_hu)';

        k_s = numel(y_s);

        fprintf('\n============================================================\n');
        fprintf('Sensitivity analysis: Excluding Hu(2021)\n');
        fprintf('Original k = %d\n',numel(idx));
        fprintf('Without Hu k = %d\n',k_s);
        fprintf('============================================================\n');

        %% -----------------------------------------------------------
        % DL
        % -----------------------------------------------------------

        w_fixed_s = 1 ./ v_s;

        theta_fixed_s = ...
            sum(w_fixed_s .* y_s) / ...
            sum(w_fixed_s);

        Q_s = ...
            sum(w_fixed_s .* ...
            (y_s - theta_fixed_s).^2);

        df_s = k_s - 1;

        C_s = ...
            sum(w_fixed_s) - ...
            sum(w_fixed_s.^2) / ...
            sum(w_fixed_s);

        if k_s > 1 && C_s > 0

            tau2_DL_s = max( ...
                0, ...
                (Q_s - df_s) / C_s);

        else

            tau2_DL_s = 0;

        end

        [theta_DL_s,var_theta_DL_s] = ...
            pooled_given_tau( ...
            y_s,v_s,tau2_DL_s);

        se_DL_s = sqrt(var_theta_DL_s);

        ci_DL_s = ...
            theta_DL_s + ...
            norminv([0.025 0.975]) .* se_DL_s;

        if Q_s > 0 && k_s > 1

            I2_s = max( ...
                0, ...
                ((Q_s-df_s)/Q_s)*100);

        else

            I2_s = 0;

        end
% ============================================================
% REML for sensitivity analysis
% ============================================================

        if k_s >= 1

            opt = optimset( ...
                'TolX', 1e-10, ...
                'MaxIter', 1000, ...
                'Display', 'off');

            upper_bound_s = max([ ...
                100 * max(v_s), ...
                100 * var(y_s), ...
                1]);

            fun_s = @(t) restrictedLogLik(t, y_s, v_s);

            try

                tau2_REML_s = fminbnd( ...
                    fun_s, ...
                    0, ...
                    upper_bound_s, ...
                    opt);

                tau2_REML_s = max(0, tau2_REML_s);

            catch ME

                warning( ...
                    'REML sensitivity optimization failed: %s. tau2_REML_s set to 0.', ...
                    E.message);

                tau2_REML_s = 0;

            end

            [theta_REML_s, var_theta_REML_s] = ...
                pooled_given_tau(y_s, v_s, tau2_REML_s);

            se_theta_REML_s = sqrt(var_theta_REML_s);

            ci_REML_s = theta_REML_s + ...
                norminv([0.025 0.975]) * se_theta_REML_s;

        else

            tau2_REML_s = NaN;
            theta_REML_s = NaN;
        ci_REML_s = [NaN NaN];

        end
        
        [theta_REML_s,var_theta_REML_s] = ...
            pooled_given_tau( ...
            y_s,v_s,tau2_REML_s);

        se_REML_s = ...
            sqrt(var_theta_REML_s);

        ci_REML_s = ...
            theta_REML_s + ...
            norminv([0.025 0.975]) .* se_REML_s;

        %% -----------------------------------------------------------
        % Prediction intervals after excluding Hu
        % -----------------------------------------------------------

        PI_DL_s = [NaN NaN];
        PI_REML_s = [NaN NaN];

        if k_s >= 3

            if k_s < 5

                warning( ...
                    ['Sensitivity PI uses only k=%d studies. ' ...
                     'Interpret cautiously.'],k_s);

            end

            df_PI_s = k_s - 2;

            tcrit_s = ...
                tinv(0.975,df_PI_s);

            PI_DL_s = ...
                theta_DL_s + ...
                [-1 1] .* ...
                tcrit_s .* ...
                sqrt(var_theta_DL_s + tau2_DL_s);

            PI_REML_s = ...
                theta_REML_s + ...
                [-1 1] .* ...
                tcrit_s .* ...
                sqrt(var_theta_REML_s + tau2_REML_s);

        end

        %% -----------------------------------------------------------
        % Sensitivity table
        % -----------------------------------------------------------

        sens_rows = { ...
            'Original_with_Hu_DL', ...
            theta_DL, ...
            S.ci_DL(1), ...
            S.ci_DL(2), ...
            S.PI_DL(1), ...
            S.PI_DL(2), ...
            S.tau2_DL, ...
            S.Q, ...
            S.I2;

            'Original_with_Hu_REML', ...
            theta_REML, ...
            S.ci_REML(1), ...
            S.ci_REML(2), ...
            S.PI_REML(1), ...
            S.PI_REML(2), ...
            S.tau2_REML, ...
            NaN, ...
            NaN;

            'Exclude_Hu_DL', ...
            theta_DL_s, ...
            ci_DL_s(1), ...
            ci_DL_s(2), ...
            PI_DL_s(1), ...
            PI_DL_s(2), ...
            tau2_DL_s, ...
            Q_s, ...
            I2_s;

            'Exclude_Hu_REML', ...
            theta_REML_s, ...
            ci_REML_s(1), ...
            ci_REML_s(2), ...
            PI_REML_s(1), ...
            PI_REML_s(2), ...
            tau2_REML_s, ...
            NaN, ...
            NaN};

        T_sens = cell2table( ...
            sens_rows, ...
            'VariableNames',{ ...
            'Scenario', ...
            'Theta', ...
            'CI_lo', ...
            'CI_hi', ...
            'PI_lo', ...
            'PI_hi', ...
            'tau2', ...
            'Q', ...
            'I2'});

        writetable( ...
            T_sens, ...
            fullfile( ...
            outdir, ...
            'sensitivity_exclude_Hu2021.xlsx'));

        %% -----------------------------------------------------------
        % Percentage change
        % -----------------------------------------------------------

        PercentChange = table( ...
            {'Theta_DL'; ...
             'tau2_DL'; ...
             'I2_DL'; ...
             'Theta_REML'; ...
             'tau2_REML'}, ...
            [ ...
            100*(theta_DL_s-S.theta_DL)/S.theta_DL; ...
            100*(tau2_DL_s-S.tau2_DL)/S.tau2_DL; ...
            100*(I2_s-S.I2)/S.I2; ...
            100*(theta_REML_s-S.theta_REML)/S.theta_REML; ...
            100*(tau2_REML_s-S.tau2_REML)/S.tau2_REML], ...
            'VariableNames', ...
            {'Metric','PercentChange'});

        writetable( ...
            PercentChange, ...
            fullfile( ...
            outdir, ...
            'sensitivity_percent_change_Hu2021.xlsx'));

        %% -----------------------------------------------------------
        % Console
        % -----------------------------------------------------------

        fprintf('\nSensitivity results:\n');

        fprintf( ...
            'DL original : %.3f [%.3f, %.3f]\n', ...
            S.theta_DL,S.ci_DL(1),S.ci_DL(2));

        fprintf( ...
            'DL no Hu    : %.3f [%.3f, %.3f]\n', ...
            theta_DL_s,ci_DL_s(1),ci_DL_s(2));

        fprintf( ...
            'REML original : %.3f [%.3f, %.3f]\n', ...
            S.theta_REML,S.ci_REML(1),S.ci_REML(2));

        fprintf( ...
            'REML no Hu    : %.3f [%.3f, %.3f]\n', ...
            theta_REML_s,ci_REML_s(1),ci_REML_s(2));

        fprintf('\n');

        fprintf( ...
            'DL tau2: %.6f -> %.6f\n', ...
            S.tau2_DL,tau2_DL_s);

        fprintf( ...
            'REML tau2: %.6f -> %.6f\n', ...
            S.tau2_REML,tau2_REML_s);

        %% -----------------------------------------------------------
        % Forest plots without Hu
        % -----------------------------------------------------------

        S_s = struct();

        S_s.group = 'Elderly_noHu';
        S_s.k = k_s;

        S_s.studies = studies(idx_no_hu);

        S_s.y = y_s(:);
        S_s.v = v_s(:);

        S_s.theta_DL = theta_DL_s;
        S_s.ci_DL = ci_DL_s;
        S_s.tau2_DL = tau2_DL_s;

        S_s.theta_REML = theta_REML_s;
        S_s.ci_REML = ci_REML_s;
        S_s.tau2_REML = tau2_REML_s;

        S_s.PI_DL = PI_DL_s;
        S_s.PI_REML = PI_REML_s;

        S_s.weights_pct_DL = ...
            100*(1./(v_s+tau2_DL_s)) ./ ...
            sum(1./(v_s+tau2_DL_s));

        S_s.weights_pct_REML = ...
            100*(1./(v_s+tau2_REML_s)) ./ ...
            sum(1./(v_s+tau2_REML_s));

        makeForestGroup( ...
            'DL', ...
            S_s, ...
            outdir);

        makeForestGroup( ...
            'REML', ...
            S_s, ...
            outdir);

    end
end

%% ========================================================================
% FINISH
% ========================================================================

fprintf('\n============================================================\n');
fprintf('META-ANALYSIS COMPLETED\n');
fprintf('Results saved in:\n%s\n',outdir);
fprintf('============================================================\n');


%% ========================================================================
% LOCAL FUNCTION 1. ajustada por el tipo de analisis
% REML RESTRICTED LOG-LIKELIHOOD
% ========================================================================

function nll = restrictedLogLik(tau2, y, v)
% ============================================================
% Restricted log-likelihood for REML estimation of tau^2
% Random-effects meta-analysis with a single pooled intercept
%
% y   = study-level effect estimates
% v   = within-study sampling variances
% tau2 = between-study variance
%
% The function returns the NEGATIVE restricted log-likelihood,
% because fminbnd performs minimization.
% ============================================================

    % Ensure column vectors
    y = y(:);
    v = v(:);

    % Invalid values
    if tau2 < 0 || any(v <= 0) || any(~isfinite(v)) || any(~isfinite(y))
        nll = Inf;
        return;
    end

    % Total variance for each study
    sigma2 = v + tau2;

    % Random-effects weights
    w = 1 ./ sigma2;

    % Sum of weights
    W = sum(w);

    % Estimated pooled effect for this value of tau2
    mu_hat = sum(w .* y) / W;

    % Weighted residual sum of squares
    residuals = y - mu_hat;
    Q_tau = sum(w .* residuals.^2);

    % Restricted negative log-likelihood
    %
    % REML for a model with one fixed effect (the pooled mean):
    %
    % -2 log L_REML =
    %     sum(log(sigma2))
    %     + log(W)
    %     + Q_tau
    %
    % Constants independent of tau2 are omitted.
    
    nll = 0.5 * ( ...
        sum(log(sigma2)) + ...
        log(W) + ...
        Q_tau );

    % Numerical protection
    if ~isfinite(nll)
        nll = Inf;
    end
end
% function out = restrictedLogLik_REML(tau2,y,v)
% 
%     % Ensure non-negative tau2
%     tau2 = max(tau2,0);
% 
%     % Total variance
%     vi = v + tau2;
% 
%     if any(~isfinite(vi)) || any(vi <= 0)
% 
%         out = Inf;
%         return;
% 
%     end
% 
%     % Random-effects weights
%     w = 1 ./ vi;
% 
%     % Sum of weights
%     W = sum(w);
% 
%     if ~isfinite(W) || W <= 0
% 
%         out = Inf;
%         return;
% 
%     end
% 
%     % Weighted pooled estimate
%     theta = ...
%         sum(w .* y) / W;
% 
%     % Weighted residual heterogeneity
%     Q_tau = ...
%         sum(w .* (y-theta).^2);
% 
%     k = numel(y);
% 
%     % Need at least 2 studies
%     if k < 2 || Q_tau <= 0
% 
%         out = Inf;
%         return;
% 
%     end
% 
%     % REML objective function
%     %
%     % Constants independent of tau2 are omitted.
%     %
%     % -2 log L_REML =
%     % sum(log(vi))
%     % + log(sum(w))
%     % + (k-1)*log(Q_tau)
% 
%     out = 0.5 * ( ...
%         sum(log(vi)) + ...
%         log(W) + ...
%         (k-1)*log(Q_tau) );
% 
%     if ~isfinite(out)
% 
%         out = Inf;
% 
%     end
% 
% end
% 
% 
%% ========================================================================
% LOCAL FUNCTION 2
% POOLED EFFECT GIVEN TAU2
% ========================================================================

function [theta,var_theta] = ...
    pooled_given_tau(y,v,tau2)

    w = 1 ./ (v + tau2);

    theta = ...
        sum(w .* y) / sum(w);

    var_theta = ...
        1 / sum(w);

end


%% ========================================================================
% LOCAL FUNCTION 3
% FOREST PLOT
% ========================================================================

function makeForestGroup(method,G,outdir)

    studies = G.studies;
    y = G.y(:);
    v = G.v(:);

    k = G.k;

    %% ---------------------------------------------------------------
    % Select estimator
    % ---------------------------------------------------------------

    if strcmpi(method,'DL')

        theta = G.theta_DL;
        ci = G.ci_DL;
        tau2 = G.tau2_DL;

        if isfield(G,'weights_pct_DL')

            weights_pct = ...
                G.weights_pct_DL;

        else

            weights_pct = ...
                100*(1./(v+tau2)) ./ ...
                sum(1./(v+tau2));

        end

    else

        theta = G.theta_REML;
        ci = G.ci_REML;
        tau2 = G.tau2_REML;

        if isfield(G,'weights_pct_REML')

            weights_pct = ...
                G.weights_pct_REML;

        else

            weights_pct = ...
                100*(1./(v+tau2)) ./ ...
                sum(1./(v+tau2));

        end

    end

    %% ---------------------------------------------------------------
    % Figure
    % ---------------------------------------------------------------

    fig = figure( ...
        'Visible','off', ...
        'Units','pixels', ...
        'Position',[100 100 1100 650]);

    ax = axes('Parent',fig);

    hold(ax,'on');

    %% ---------------------------------------------------------------
    % Plot study effects
    % ---------------------------------------------------------------

    ypos = (1:k)';

    max_w = max(weights_pct);

    if max_w <= 0
        max_w = 1;
    end

    for i = 1:k

        xi = y(i);

        sei = sqrt(v(i));

        ci_low = ...
            xi - 1.96*sei;

        ci_high = ...
            xi + 1.96*sei;

        % Confidence interval
        plot( ...
            ax, ...
            [ci_low ci_high], ...
            [ypos(i) ypos(i)], ...
            'k-', ...
            'LineWidth',1.4);

        % Study marker
        sz = ...
            max(25, ...
            150*weights_pct(i)/max_w);

        scatter( ...
            ax, ...
            xi, ...
            ypos(i), ...
            sz, ...
            'k', ...
            'filled');

    end

    %% ---------------------------------------------------------------
    % X limits
    % ---------------------------------------------------------------

    xleft = ...
        min([y - 3*sqrt(v); ci(1)]);

    xright = ...
        max([y + 3*sqrt(v); ci(2)]);

    if ~isfinite(xleft) || ...
       ~isfinite(xright) || ...
       xleft == xright

        xleft = min(y)-1;
        xright = max(y)+1;

    end

    xmargin = ...
        0.12*(xright-xleft);

    if xmargin == 0

        xmargin = ...
            max(1,0.1*abs(xright));

    end

    xlim(ax, ...
        [xleft-xmargin xright+xmargin]);

    %% ---------------------------------------------------------------
    % Null line
    % ---------------------------------------------------------------

    yl = [0 k+1];

    plot( ...
        ax, ...
        [0 0], ...
        yl, ...
        '--', ...
        'Color',[0.5 0.5 0.5], ...
        'LineWidth',1);

    %% ---------------------------------------------------------------
    % Pooled diamond
    % ---------------------------------------------------------------

    diamond_y = 0.35;

    diamond_h = 0.45;

    patch_x = [ ...
        ci(1), ...
        theta, ...
        ci(2), ...
        theta];

    patch_y = ...
        diamond_y + ...
        [-diamond_h/2, ...
          0, ...
          diamond_h/2, ...
          0];

    patch( ...
        ax, ...
        patch_x, ...
        patch_y, ...
        'k', ...
        'FaceAlpha',0.4, ...
        'EdgeColor','k');

    %% ---------------------------------------------------------------
    % Y-axis
    % ---------------------------------------------------------------

    set( ...
        ax, ...
        'YTick',ypos, ...
        'YTickLabel',studies, ...
        'YDir','reverse');

    %% ---------------------------------------------------------------
    % Pooled annotation
    % ---------------------------------------------------------------

    txtPooled = sprintf( ...
        'Pooled = %.2f [%.2f, %.2f]', ...
        theta,ci(1),ci(2));

    text( ...
        ax, ...
        mean(ci), ...
        diamond_y+0.55, ...
        txtPooled, ...
        'HorizontalAlignment','center', ...
        'FontWeight','bold');

    %% ---------------------------------------------------------------
    % Study estimates and weights
    % ---------------------------------------------------------------

    xlims = xlim(ax);

    xtext = ...
        xlims(2) - ...
        0.015*range(xlims);

    for i = 1:k

        txt = sprintf( ...
            '%.2f [%.2f, %.2f]  (%.1f%%)', ...
            y(i), ...
            y(i)-1.96*sqrt(v(i)), ...
            y(i)+1.96*sqrt(v(i)), ...
            weights_pct(i));

        text( ...
            ax, ...
            xtext, ...
            ypos(i), ...
            txt, ...
            'HorizontalAlignment','right', ...
            'FontSize',9);

    end

    %% ---------------------------------------------------------------
    % Labels
    % ---------------------------------------------------------------

    xlabel( ...
        ax, ...
        'Mean Difference (Control - Sarcopenia)');

    if isfield(G,'group')

        groupTitle = G.group;

    else

        groupTitle = 'Meta-analysis';

    end

    title( ...
        ax, ...
        sprintf('%s - %s', ...
        groupTitle,upper(method)), ...
        'Interpreter','none');

    grid(ax,'on');

    %% ---------------------------------------------------------------
    % Save
    % ---------------------------------------------------------------

    safeGroup = ...
        matlab.lang.makeValidName(groupTitle);

    fname = fullfile( ...
        outdir, ...
        sprintf( ...
        'Forest_%s_%s', ...
        safeGroup, ...
        upper(method)));

    try

        exportgraphics( ...
            fig, ...
            [fname '.png'], ...
            'Resolution',300, ...
            'ContentType','image');

        exportgraphics( ...
            fig, ...
            [fname '.pdf'], ...
            'ContentType','vector');

    catch

        saveas( ...
            fig, ...
            [fname '.png']);

        saveas( ...
            fig, ...
            [fname '.pdf']);

    end

    close(fig);

end