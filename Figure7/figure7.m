clear; clc; close all;

% ---------------------- Paths & options ----------------------
% Update this path if needed:
csv_file  = "markov_clinical.csv";

dt     = 1;
p0_3   = [0.25; 0.25];  % initial for early overlays
opts   = odeset('RelTol',1e-8,'AbsTol',1e-10');
FIG_W_IN = 6.0;
FIG_H_IN = 4.5;
AX_FONT = 16;  % larger axes fonts [[memory:6659181]]

% ---------------------- Load clinical matrix ----------------------
C        = readcell(csv_file, 'Delimiter', ',');  
rowNames = string(C(:,1));
vals     = C(:,2:end);                           

get_row_cells = @(name) vals(find(strcmpi(rowNames, name), 1, 'first'), :);
to_num_vec    = @(cells) cellfun(@(x) local_to_double(x), cells);

id_vec      = to_num_vec(get_row_cells("ID"));
relapse_vec = to_num_vec(get_row_cells("Relapse"));
dstate_vec  = string(get_row_cells("DiseaseState"));

uIDs = unique(id_vec(~isnan(id_vec)));

% ---------------------- Per-ID params via lumping ----------------------
ID_list    = [];
alpha_list = []; beta_list = []; gamma_list = []; delta_list = []; sigma_list = []; tau_list = [];
grp_relapse= false(0,1);
grp_remiss = false(0,1);

for k = 1:numel(uIDs)
    IDk     = uIDs(k);
    col_idx = find(id_vec == IDk);
    if isempty(col_idx), continue; end

    % Build mean 4x4 P per ID
    M  = zeros(4,4);
    okM= true;
    for i = 1:4
        for j = 1:4
            rn = "M"+i+""+j;
            row_cells = get_row_cells(rn);
            if isempty(row_cells), okM=false; break; end
            row_vals  = to_num_vec(row_cells);
            M(i,j)    = mean(row_vals(col_idx), 'omitnan');
        end
        if ~okM, break; end
    end
    if ~okM, continue; end

    P = M; 
    Q    = to_generator(P, dt);
    pars = lump_to_three(Q, P);

    ID_list(end+1,1)    = IDk;
    alpha_list(end+1,1) = pars.alpha;
    beta_list(end+1,1)  = pars.beta;
    gamma_list(end+1,1) = pars.gamma;
    delta_list(end+1,1) = pars.delta;
    sigma_list(end+1,1) = pars.sigma;
    tau_list(end+1,1)   = pars.tau;

    this_rel = relapse_vec(col_idx);
    grp_relapse(end+1,1) = any(this_rel == 1);

    this_ds  = dstate_vec(col_idx);
    grp_remiss(end+1,1)  = any(strcmpi(this_ds, "Remission"));
end

% ---------------------- Colors ----------------------
relapseColor    = hex2rgb('#86ACD8');   
remissionColor  = hex2rgb('#D27BDD');  
diffColor       = hex2rgb('#77DD77');      % differentiation promoter (σ↑, β↑)
dediffColor     = hex2rgb('#FF6961');      % de-differentiation inhibitor (α↓, τ↓)

% ---------------------- Overlayed time series (group means) ----------------------
t_eval = linspace(0,10,201);
nT  = numel(t_eval);
nID = numel(ID_list);
P1 = NaN(nT, nID); P3 = NaN(nT, nID);

for r = 1:numel(ID_list)
    pars.alpha = alpha_list(r); pars.beta  = beta_list(r);
    pars.gamma = gamma_list(r); pars.delta = delta_list(r);
    pars.sigma = sigma_list(r); pars.tau   = tau_list(r);
    [~, y]     = ode45(@(t,y) rhs_3state_local(t,y,pars), t_eval, p0_3, opts);
    P1(:,r) = y(:,1); P3(:,r) = y(:,2);
end

col_p1_base = hex2rgb('#A976BB'); col_p3_base = hex2rgb('#F6AFCC');
col_p1_bright = lightenColor(col_p1_base, 0.5);
col_p3_bright = lightenColor(col_p3_base, 0.5);
col_p1_mean   = darkenColor(col_p1_base, 0.2);
col_p3_mean   = darkenColor(col_p3_base, 0.2);

idxR = find(grp_remiss);
if ~isempty(idxR)
    fRem = figure('Color','w'); hold on;box on;
    for k = idxR(:)'
        plot(t_eval, P1(:,k), 'Color', col_p1_bright, 'LineWidth', 0.8);
        plot(t_eval, P3(:,k), 'Color', col_p3_bright, 'LineWidth', 0.8);
    end
    p1m = mean(P1(:,idxR), 2, 'omitnan'); p3m = mean(P3(:,idxR), 2, 'omitnan');
    h1 = plot(t_eval, p1m, 'Color', col_p1_mean, 'LineWidth', 3, 'DisplayName','CD34+/CD38-');
    h3 = plot(t_eval, p3m, 'Color', col_p3_mean, 'LineWidth', 3, 'DisplayName','CD34-/CD38+');
    grid on; ylim([0 1]); xlim([0 2]); legend([h1 h3],'Location','best');
    ax=gca; ax.FontSize=AX_FONT; ax.FontWeight='bold';
    set_figure_size(fRem, FIG_W_IN, FIG_H_IN);
    exportgraphics(fRem,'remission_overlay.png','Resolution',300);
end

idxL = find(grp_relapse);
if ~isempty(idxL)
    fRel = figure('Color','w'); hold on;box on;
    for k = idxL(:)'
        plot(t_eval, P1(:,k), 'Color', col_p1_bright, 'LineWidth', 0.8);
        plot(t_eval, P3(:,k), 'Color', col_p3_bright, 'LineWidth', 0.8);
    end
    p1m = mean(P1(:,idxL), 2, 'omitnan'); p3m = mean(P3(:,idxL), 2, 'omitnan');
    h1 = plot(t_eval, p1m, 'Color', col_p1_mean, 'LineWidth', 3, 'DisplayName','CD34+/CD38-');
    h3 = plot(t_eval, p3m, 'Color', col_p3_mean, 'LineWidth', 3, 'DisplayName','CD34-/CD38+');
    grid on; ylim([0 1]); xlim([0 2]); legend([h1 h3],'Location','best');
    ax=gca; ax.FontSize=AX_FONT; ax.FontWeight='bold';
    set_figure_size(fRel, FIG_W_IN, FIG_H_IN);
    exportgraphics(fRel,'relapse_overlay.png','Resolution',300);
end

% ---------------------- Param differences: Relapse vs Remission ----------------------
d_sigma_tau = sigma_list + beta_list - tau_list - alpha_list;   % (σ + β) - (τ + α)

idxRel = idxL; idxRem = idxR;

fDiff = figure('Color','w');
tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
nexttile; hold on;box on;
valsRel = d_sigma_tau(idxRel);
valsRem = d_sigma_tau(idxRem);
x1 = categorical(repmat("Relapse",  numel(valsRel),1), ["Relapse","Remission"]);
x2 = categorical(repmat("Remission",numel(valsRem),1), ["Relapse","Remission"]);
if ~isempty(valsRel), boxchart(x1, valsRel, 'BoxFaceAlpha',0.6, 'BoxFaceColor', relapseColor,   'LineWidth',1.2); end
if ~isempty(valsRem), boxchart(x2, valsRem, 'BoxFaceAlpha',0.6, 'BoxFaceColor', remissionColor, 'LineWidth',1.2); end
yline(0,'k:','LineWidth',1);
ax=gca; ax.FontSize=AX_FONT; ax.FontWeight='bold';
set_figure_size(fDiff, FIG_W_IN, FIG_H_IN);
exportgraphics(fDiff,'param_differences_relapse_vs_remission.png','Resolution',300);

% ---------------------- Remission composition targets ----------------------
p1_target = 0.06;   % CD34+/CD38- in BM Remission (Fig. 1H)
p3_target = 0.58;   % CD34-/CD38+ in BM Remission (Fig. 1H)
p0_3 = [p1_target; p3_target];

% ---------------------- Intervention factors ----------------------
alpha_tau_factor = 1.0;  % used for combo
sigma_factor = 2.5;      
beta_factor  = 2.5;      

diff_sigma_factor = 2.5;   
diff_beta_factor  = 2.5;   
dediff_factor     = 0.1;   

% ---------------------- Endpoints for phase/box figures (relapse group) ----------------------
P1_end_base   = nan(numel(idxL),1); P3_end_base   = nan(numel(idxL),1);
P1_end_diff   = nan(numel(idxL),1); P3_end_diff   = nan(numel(idxL),1);
P1_end_dediff = nan(numel(idxL),1); P3_end_dediff = nan(numel(idxL),1);

for ii = 1:numel(idxL)
    r = idxL(ii);
    base.alpha = alpha_list(r); base.beta  = beta_list(r);
    base.gamma = gamma_list(r); base.delta = delta_list(r);
    base.sigma = sigma_list(r); base.tau   = tau_list(r);

	[~, y] = ode45(@(t,y) rhs_3state_local(t,y,base), t_eval, p0_3, opts);
	P1_end_base(ii) = y(end,1); P3_end_base(ii) = y(end,2);

    % Differentiation promoter = σ↑ + β↑
    pDiff = base; 
    pDiff.sigma = diff_sigma_factor*base.sigma;
    pDiff.beta  = diff_beta_factor *base.beta;
    [~, yDiff] = ode45(@(t,y) rhs_3state_local(t,y,pDiff), t_eval, p0_3, opts);
    P1_end_diff(ii) = yDiff(end,1); P3_end_diff(ii) = yDiff(end,2);

    % De-differentiation inhibitor = α↓ + τ↓
    pDe = base;
    pDe.alpha = dediff_factor*base.alpha;
    pDe.tau   = dediff_factor*base.tau;
    [~, yDe] = ode45(@(t,y) rhs_3state_local(t,y,pDe), t_eval, p0_3, opts);
    P1_end_dediff(ii) = yDe(end,1); P3_end_dediff(ii) = yDe(end,2);
end

% ---------------------- Relapse (combo σ↑+β↑) overlay (single plot) ----------------------
if ~isempty(idxL)
    p0_init = [0.25; 0.25];
    P1L = NaN(nT, numel(idxL)); P3L = NaN(nT, numel(idxL));
    for j = 1:numel(idxL)
        r = idxL(j);
        base.alpha = alpha_list(r); base.beta  = beta_list(r);
        base.gamma = gamma_list(r); base.delta = delta_list(r);
        base.sigma = sigma_list(r); base.tau   = tau_list(r);
        pC = base; pC.sigma = sigma_factor*base.sigma; pC.beta = beta_factor*base.beta; 
        pC.alpha = alpha_tau_factor * base.alpha;
        pC.tau   = alpha_tau_factor * base.tau;
        [~, yC] = ode45(@(t,y) rhs_3state_local(t,y,pC), t_eval, p0_init, opts);
        P1L(:,j) = yC(:,1); P3L(:,j) = yC(:,2);
    end
    p1m_L = mean(P1L,2,'omitnan'); p3m_L = mean(P3L,2,'omitnan');

    fRelSingle = figure('Color','w'); hold on; box on;
    for j = 1:size(P1L,2), plot(t_eval, P1L(:,j), 'Color', col_p1_bright, 'LineWidth', 0.8); end
    for j = 1:size(P3L,2), plot(t_eval, P3L(:,j), 'Color', col_p3_bright, 'LineWidth', 0.8); end
    h1 = plot(t_eval, p1m_L, 'Color', col_p1_mean, 'LineWidth', 3, 'DisplayName','CD34+/CD38-');
    h3 = plot(t_eval, p3m_L, 'Color', col_p3_mean, 'LineWidth', 3, 'DisplayName','CD34-/CD38+');

    ylim([0 1]); xlim([0 2]); grid on;
    legend([h1 h3],'Location','best');
    ax=gca; ax.FontSize=AX_FONT; ax.FontWeight='bold';
    set_figure_size(fRelSingle, FIG_W_IN, FIG_H_IN);
    exportgraphics(fRelSingle,'relapse_overlay_combo_singleplot.png','Resolution',300);
end

% ---------------------- Phase plot: baseline → interventions (relapse group) ----------------------
if ~isempty(idxL)
	comboColor = [0.20 0.20 0.20];
	baseEdge   = darkenColor(relapseColor, 0.50);

	fPhase = figure('Color','w'); hold on; box on;
	s1 = scatter(P1_end_base,  P3_end_base, 30, relapseColor, 'filled', 'DisplayName','Baseline (relapse)');
	s1.MarkerEdgeColor = baseEdge;

	% Endpoints for two interventions
	s2 = scatter(P1_end_diff,    P3_end_diff,    40, diffColor,   'filled', 'DisplayName','Diff promoter (σ↑, β↑)');
	s3 = scatter(P1_end_dediff,  P3_end_dediff,  40, dediffColor, 'filled', 'DisplayName','De-diff inhibitor (α↓, τ↓)');

	% Arrows from baseline → each intervention
	for ii = 1:numel(P1_end_base)
		if all(isfinite([P1_end_base(ii) P3_end_base(ii) P1_end_diff(ii) P3_end_diff(ii)]))
			quiver(P1_end_base(ii), P3_end_base(ii), ...
			       P1_end_diff(ii)-P1_end_base(ii), P3_end_diff(ii)-P3_end_base(ii), ...
			       0, 'Color', diffColor, 'LineWidth',0.9, 'MaxHeadSize',0.7);
		end
		if all(isfinite([P1_end_base(ii) P3_end_base(ii) P1_end_dediff(ii) P3_end_dediff(ii)]))
			quiver(P1_end_base(ii), P3_end_base(ii), ...
			       P1_end_dediff(ii)-P1_end_base(ii), P3_end_dediff(ii)-P3_end_base(ii), ...
			       0, 'Color', dediffColor, 'LineWidth',0.9, 'MaxHeadSize',0.7);
		end
	end

	plot(p1_target, p3_target, 'p', 'MarkerSize', 12, ...
	     'MarkerFaceColor', remissionColor, 'MarkerEdgeColor', 'k', ...
	     'DisplayName','Remission target');
	xline(p1_target,'k:','LineWidth',2);
	yline(p3_target,'k:','LineWidth',2);

    axis([0 1 0 1]);  grid on;
    ax=gca; ax.FontSize=AX_FONT; ax.FontWeight='bold';
    set_figure_size(fPhase, FIG_W_IN, FIG_H_IN);
    exportgraphics(fPhase,'relapse_to_remission_phase.png','Resolution',300);
end

% ---------------------- p1 boxplots: Remission vs Relapse vs interventions ----------------------
p1_end_rem_base = NaN(numel(idxR),1);
if ~isempty(idxR)
	for ii = 1:numel(idxR)
		r = idxR(ii);
		base.alpha = alpha_list(r); base.beta  = beta_list(r);
		base.gamma = gamma_list(r); base.delta = delta_list(r);
		base.sigma = sigma_list(r); base.tau   = tau_list(r);
		[~, y] = ode45(@(t,y) rhs_3state_local(t,y,base), t_eval, p0_3, opts);
		p1_end_rem_base(ii) = y(end,1);
	end
end

p1_rem_base = p1_end_rem_base(~isnan(p1_end_rem_base));
p1_rel_base2 = P1_end_base(~isnan(P1_end_base));
p1_rel_diff  = P1_end_diff(~isnan(P1_end_diff));
p1_rel_ded   = P1_end_dediff(~isnan(P1_end_dediff));

fP1four = figure('Color','w'); hold on; box on;
cats = ["Remission","Relapse","Diff promoter","De-diff inhibitor"];
if ~isempty(p1_rem_base),  boxchart(categorical(repmat("Remission",numel(p1_rem_base),1), cats),   p1_rem_base,  'BoxFaceAlpha',0.6, 'BoxFaceColor', remissionColor, 'LineWidth',1.2); end
if ~isempty(p1_rel_base2), boxchart(categorical(repmat("Relapse",  numel(p1_rel_base2),1), cats), p1_rel_base2, 'BoxFaceAlpha',0.6, 'BoxFaceColor', relapseColor,   'LineWidth',1.2); end
if ~isempty(p1_rel_diff),  boxchart(categorical(repmat("Diff promoter",numel(p1_rel_diff),1), cats),  p1_rel_diff,  'BoxFaceAlpha',0.8, 'BoxFaceColor', diffColor,     'LineWidth',1.2); end
if ~isempty(p1_rel_ded),   boxchart(categorical(repmat("De-diff inhibitor",numel(p1_rel_ded),1), cats), p1_rel_ded, 'BoxFaceAlpha',0.8, 'BoxFaceColor', dediffColor,   'LineWidth',1.2); end
ylim([0 1]); grid on;
ax = gca;
ax.XTickLabel = [];
ax.FontSize = AX_FONT;
ax.FontWeight = 'bold';
set_figure_size(fP1four, FIG_W_IN, FIG_H_IN);
exportgraphics(fP1four,'p1_relapse_vs_remission_combo.png','Resolution',300);

% ---------------------- Relapse (De-diff inhibitor α↓+τ↓) overlay (single plot) ----------------------
if ~isempty(idxL)
	p0_init = [0.25; 0.25];
	P1L = NaN(nT, numel(idxL)); P3L = NaN(nT, numel(idxL));
	for j = 1:numel(idxL)
		r = idxL(j);
		base.alpha = alpha_list(r); base.beta  = beta_list(r);
		base.gamma = gamma_list(r); base.delta = delta_list(r);
		base.sigma = sigma_list(r); base.tau   = tau_list(r);

		pDe = base;
		pDe.alpha = dediff_factor*base.alpha;
		pDe.tau   = dediff_factor*base.tau;
		[~, yE] = ode45(@(t,y) rhs_3state_local(t,y,pDe), t_eval, p0_init, opts);
		P1L(:,j) = yE(:,1); P3L(:,j) = yE(:,2);
	end
	p1m_L = mean(P1L,2,'omitnan'); p3m_L = mean(P3L,2,'omitnan');

	fDeSingle = figure('Color','w'); hold on; box on;
	for j = 1:size(P1L,2), plot(t_eval, P1L(:,j), 'Color', col_p1_bright, 'LineWidth', 0.8); end
	for j = 1:size(P3L,2), plot(t_eval, P3L(:,j), 'Color', col_p3_bright, 'LineWidth', 0.8); end
	h1 = plot(t_eval, p1m_L, 'Color', col_p1_mean, 'LineWidth', 3, 'DisplayName','CD34+/CD38-');
	h3 = plot(t_eval, p3m_L, 'Color', col_p3_mean, 'LineWidth', 3, 'DisplayName','CD34-/CD38+');

	ylim([0 1]); xlim([0 2]); grid on;
	legend([h1 h3],'Location','best');
	ax=gca; ax.FontSize=AX_FONT; ax.FontWeight='bold';
    set_figure_size(fDeSingle, FIG_W_IN, FIG_H_IN);
	exportgraphics(fDeSingle,'relapse_overlay_dediff_singleplot.png','Resolution',300);
end

% ============================ Helpers ============================
function v = local_to_double(x)
    if isnumeric(x)
        v = x;
    elseif isstring(x) || ischar(x)
        xx = strtrim(string(x));
        if xx == "" || strcmpi(xx,"NaN")
            v = NaN;
        else
            t = str2double(xx);
            if isnan(t), v = NaN; else, v = t; end
        end
    else
        v = NaN;
    end
end

function Q = to_generator(P, dt)
    A = logm(P');              % matrix log of P^T
    Q = real(A) / dt;
    for j=1:4
        for i=1:4
            if i~=j && Q(i,j)<0, Q(i,j)=0; end
        end
        Q(j,j) = -sum(Q([1:j-1 j+1:4], j));
    end
end

function pars = lump_to_three(Q, P)
    [V,D] = eig(P');
    [~,ix]= min(abs(diag(D)-1));
    pi = abs(V(:,ix)); pi = pi/sum(pi);
    w2 = pi(2)/(pi(2)+pi(4)); w4 = 1-w2;

    pars.alpha = w2*Q(1,2) + w4*Q(1,4); % O->1
    pars.gamma = w2*Q(3,2) + w4*Q(3,4); % O->3
    pars.beta  = Q(2,1) + Q(4,1);       % 1->O
    pars.delta = Q(2,3) + Q(4,3);       % 3->O
    pars.sigma = Q(3,1);                % 1->3
    pars.tau   = Q(1,3);                % 3->1
end

function rgb = hex2rgb(hex)
    hex = char(regexprep(hex,'^#',''));
    assert(numel(hex)==6, 'hex must be #RRGGBB');
    rgb = [hex2dec(hex(1:2)) hex2dec(hex(3:4)) hex2dec(hex(5:6))]/255;
end

function c2 = lightenColor(c, f)
    c2 = c*(1-f) + f*[1 1 1];
end

function c2 = darkenColor(c, f)
    c2 = c*(1-f);
end

function dydt = rhs_3state_local(~, y, pars)
    p1 = y(1); p3 = y(2); pO = 1 - p1 - p3;
    alpha = pars.alpha; beta = pars.beta; gamma = pars.gamma; delta = pars.delta;
    sigma = pars.sigma; tau   = pars.tau;
    dp1 = alpha*pO - beta*p1 - sigma*p1 + tau*p3;
    dp3 = gamma*pO - delta*p3 + sigma*p1 - tau*p3;
    dydt = [dp1; dp3];
end

function set_figure_size(fig, w, h)
	oldUnits = fig.Units;
	fig.Units = 'inches';
	pos = fig.Position;
	fig.Position = [pos(1) pos(2) w h];
	fig.Units = oldUnits;
end