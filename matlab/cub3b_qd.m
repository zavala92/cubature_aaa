function cub3b_qd()
%CUB3B_QD   Quadrature domain with a derivative node (Figure 4).
%   phi(D) with phi(zeta) = zeta + a*zeta^2, a = 0.3. Exact rule:
%
%       int int f dA = pi*((1 + 2a^2)*f(0) + a*f'(0)).
%
%   The fit has simple poles only, so f'(0) shows up as two close nodes
%   with large weights that nearly cancel. For degrees 2:2:20 and
%   f = exp(z) the script prints the error of
%     raw      the rule as computed
%     moments  the pair replaced by mu0*f(zC) + mu1*f'(zC)
%     refit    a second fit with 1/(z-zC), 1/(z-zC)^2 in the basis
%   and writes cub3b_qd.pdf and cub3b_validation.csv (one row per degree
%   with a detected pair).
cfg = cub_config();
a = 0.3;
phi  = @(zeta) zeta + a*zeta.^2;
phip = @(zeta) 1 + 2*a*zeta;
t = 2*pi*(0:cfg.M-1)'/cfg.M;
zq = phi(exp(1i*t)); omq = conj(zq);
zfun  = @(s) phi(exp(1i*s));
zpfun = @(s) 1i*exp(1i*s).*phip(exp(1i*s));
omfun = @(s) conj(phi(exp(1i*s)));
% diameter of the sample set, scales the pair separation in detect_cluster
diam = max(abs(zq - zq.'), [], 'all');
tv = 2*pi*((0:4*cfg.M-1)' + 0.5)/(4*cfg.M);
tv2 = 2*pi*((0:8*cfg.M-1)' + 0.5)/(8*cfg.M);
zv = phi(exp(1i*tv)); zv2 = phi(exp(1i*tv2));
exactArea = pi*(1+2*a^2); exactDerivative = pi*a;
fprintf('AAA implementation: %s\n', cfg.deps.aaa);
rows = [];

f1 = @(z) exp(z); f1p = @(z) exp(z);
Iref = cub_ref(f1, zfun, zpfun, omfun);
Iexact = pi*((1 + 2*a^2)*f1(0) + a*1);        % f'(0) = 1 for e^z
fprintf('QD: boundary ref vs exact identity: %.2e\n', abs(Iref - Iexact));

% test functions and their derivatives
fs  = {f1, @(z) 1./(1.5 - z), @(z) z.^3, @(z) z.^7, @(z) cos(3*z)};
fps = {f1p, @(z) 1./(1.5 - z).^2, @(z) 3*z.^2, @(z) 7*z.^6, @(z) -3*sin(3*z)};
Iex = @(f, fp) pi*((1 + 2*a^2)*f(0) + a*fp(0));

dn = []; de_raw = []; de_mom = []; de_ref = [];
fprintf('%4s %3s %10s %10s %10s %9s %9s %9s\n', 'deg', 'n', 'raw', ...
        'moments', 'refit', 'eta_sep', 'eta_canc', 'share');
for deg = 2:2:20
    [zk, wk, ~, zout] = cub_rule(zq, omq, deg, cfg.pdeg);
    if isempty(zk), continue; end
    dn(end+1) = deg; %#ok<*AGROW>
    de_raw(end+1) = abs(sum(wk.*f1(zk)) - Iexact)/abs(Iexact);
    [inC, zC, diag] = detect_cluster(zk, wk, diam);
    Imom = moment_eval(zk, wk, inC, zC, f1, f1p);
    de_mom(end+1) = abs(Imom - Iexact)/abs(Iexact);
    if any(inC)
        [c1, c2, zrest, crest, info] = confluent_refit( ...
            zq, omq, zC, zk(~inC), zout, cfg.pdeg, zv, zv2);
        Irefit = pi*(c1*f1(zC) + c2*f1p(zC) + sum(crest.*f1(zrest)));
        delta = zk(inC)-zC;
        % error of the moment formula for exp, i.e. the Taylor tail
        % delta^2/2! + delta^3/3! + ... about zC, summed term by term
        terms = delta.^2/2; tail = terms;
        for j=3:12, terms = terms.*delta/j; tail = tail+terms; end
        trunc = abs(exp(zC)*sum(wk(inC).*tail));
        truncBound = exp(real(zC)+max(abs(delta)))*sum(abs(wk(inC)).*abs(delta).^2)/2;
        assert(trunc <= truncBound*(1+1e-10));
        efAll = 0;
        for j=1:numel(fs)
            Ij=Iex(fs{j},fps{j});
            Qj=pi*(c1*fs{j}(zC)+c2*fps{j}(zC)+sum(crest.*fs{j}(zrest)));
            efAll=max(efAll,abs(Qj-Ij)/max(1,abs(Ij)));
        end
        sens0 = pi*(abs(c1)+sum(abs(crest))); sens1=pi*abs(c2);
        fprintf(['  validation deg %d: fit %.3e, Mv3200 %.3e, Mv6400 %.3e; ' ...
                 'trunc(exp) %.3e bound %.3e; five-f %.3e; value sensitivity %.12g derivative sensitivity %.12g\n'], ...
                 deg, info.fit, info.val, info.val2, trunc, truncBound, efAll, sens0, sens1);
        fprintf('  centre %.3e; coefficient defects %.3e %.3e\n', abs(zC), ...
                 abs(pi*c1-exactArea), abs(pi*c2-exactDerivative));
        rows(end+1,:)=[deg,info.fit,info.val,info.val2,trunc,truncBound,efAll,sens0,sens1];
    else
        Irefit = NaN;  % no pair at this degree
    end
    de_ref(end+1) = abs(Irefit - Iexact)/abs(Iexact);
    fprintf('%4d %3d %10.2e %10.2e %10.2e %9.1e %9.1e %9.3f\n', deg, ...
            numel(zk), de_raw(end), de_mom(end), de_ref(end), ...
            diag.sep, diag.cancel, diag.share);
    if deg == 10
        [~, ord] = sort(abs(wk), 'descend');
        pair = zk(ord(1:2));
        fprintf(['  deg 10: confluent pair at %.3e%+.3ei and %.3e%+.3ei ' ...
                 '(separation %.1e), |weights| %.2e, %.2e\n'], ...
                real(pair(1)), imag(pair(1)), real(pair(2)), imag(pair(2)), ...
                abs(pair(1) - pair(2)), abs(wk(ord(1))), abs(wk(ord(2))));
        mu0 = sum(wk(inC)); mu1 = sum(wk(inC).*(zk(inC) - zC));
        fprintf('  cluster centre |z_C| = %.2e; Lambda_n = %.3e\n', ...
                abs(zC), sum(abs(wk))/(pi*(1 + 2*a^2)));
        fprintf('  moments: |mu_0 - pi(1+2a^2)| = %.2e, |mu_1 - pi a| = %.2e\n', ...
                abs(mu0 - pi*(1 + 2*a^2)), abs(mu1 - pi*a));
        fprintf('  refit:   |pi c_1 - pi(1+2a^2)| = %.2e, |pi c_2 - pi a| = %.2e\n', ...
                abs(pi*c1 - pi*(1 + 2*a^2)), abs(pi*c2 - pi*a));
        % largest error over the five test functions
        er = 0; em = 0; ef = 0;
        for k = 1:numel(fs)
            Ik = Iex(fs{k}, fps{k});
            er = max(er, abs(sum(wk.*fs{k}(zk)) - Ik)/max(1, abs(Ik)));
            em = max(em, abs(moment_eval(zk, wk, inC, zC, fs{k}, fps{k}) - Ik)/max(1, abs(Ik)));
            Ir = pi*(c1*fs{k}(zC) + c2*fps{k}(zC) + sum(crest.*fs{k}(zrest)));
            ef = max(ef, abs(Ir - Ik)/max(1, abs(Ik)));
        end
        fprintf('  exactness over 5 integrands: raw %.1e, moments %.1e, refit %.1e\n', er, em, ef);
        pairkeep = pair; zCkeep = zC;
    end
end

FS = cfg.FS;
figure('Color','w','Position',cfg.fig.position);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile
plot(real(zq([1:end 1])), imag(zq([1:end 1])), '-', 'Color', cfg.col.curve, ...
     'LineWidth', cfg.plot.curveWidth); hold on
plot(real(pairkeep), imag(pairkeep), '.', 'Color', cfg.col.node, 'MarkerSize', 16)
plot(real(zCkeep), imag(zCkeep), 'k+', 'MarkerSize', 9, 'LineWidth', 1.1)
text(0.06, 0.09, 'confluent pair', 'FontSize', FS-2, 'Color', cfg.col.node)
axis equal; box on
xlabel('Re $z$', 'FontSize', FS); ylabel('Im $z$', 'FontSize', FS)
nexttile
semilogy(dn, de_raw, 'o', 'Color', cfg.col.node, ...
         'MarkerSize', cfg.plot.markerSize, 'LineWidth', 1); hold on
semilogy(dn, de_mom, 's', 'Color', cfg.col.violet, ...
         'MarkerSize', cfg.plot.markerSize, 'LineWidth', 1)
semilogy(dn, de_ref, '^', 'Color', cfg.col.gauss, ...
         'MarkerSize', cfg.plot.markerSize, 'LineWidth', 1)
ylim(cfg.plot.errorLimits); grid on; box on; set(gca,'YMinorGrid','off')
xlabel('AAA degree', 'FontSize', FS)
ylabel('relative error', 'FontSize', FS)
legend({'two-node rule', 'two-moment approximation', 'confluent refit'}, ...
       'FontSize', FS-3, 'Box', 'off', 'Location', 'northeast')
for ax = findall(gcf, 'Type', 'axes')'
    ax.Toolbar.Visible = 'off';
end
writematrix(rows,fullfile(cfg.here,'cub3b_validation.csv'));
exportgraphics(gcf, fullfile(cfg.here, 'cub3b_qd.pdf'), ...
               'ContentType', 'vector', 'BackgroundColor', 'white');
end

function [inC, zC, diag] = detect_cluster(zk, wk, diam)
% Pair = the two nodes of largest |weight|. It is accepted if the nodes
% are closer than 1e-3*diam and |w1 + w2| < 1e-2*(|w1| + |w2|). zC is
% their mean, weighted by |w|.
inC = false(size(zk)); zC = NaN;
diag = struct('sep', NaN, 'cancel', NaN, 'share', NaN);
if numel(zk) < 2, return; end
[~, ord] = sort(abs(wk), 'descend');
p = ord(1:2);
sep = abs(zk(p(1)) - zk(p(2)))/diam;
canc = abs(sum(wk(p)))/sum(abs(wk(p)));
share = sum(abs(wk(p)))/sum(abs(wk));
diag = struct('sep', sep, 'cancel', canc, 'share', share);
if sep < 1e-3 && canc < 1e-2
    inC(p) = true;
    zC = sum(abs(wk(p)).*zk(p))/sum(abs(wk(p)));
end
end

function I = moment_eval(zk, wk, inC, zC, f, fp)
% Rule with the pair replaced by mu0*f(zC) + mu1*f'(zC).
if ~any(inC), I = sum(wk.*f(zk)); return; end
mu0 = sum(wk(inC)); mu1 = sum(wk(inC).*(zk(inC) - zC));
I = mu0*f(zC) + mu1*fp(zC) + sum(wk(~inC).*f(zk(~inC)));
end

function [c1, c2, zrest, crest, info] = confluent_refit(zs, om, zC, zrest, zout, pdeg, zv, zv2)
% Fit again with 1/(z-zC), 1/(z-zC)^2 in place of the pair. info has the
% residual on the fitting points and on the grids zv, zv2.
sc = max(abs(zs));
A = [1./(zs - zC), 1./(zs - zC).^2, 1./(zs - zrest.'), 1./(zs - zout.'), ...
     (zs/sc).^(0:pdeg)];
c = cub_lsq_svd(A, om);
c1 = c(1); c2 = c(2); crest = c(2 + (1:numel(zrest)));
Av = [1./(zv-zC), 1./(zv-zC).^2, 1./(zv-zrest.'), 1./(zv-zout.'), (zv/sc).^(0:pdeg)];
Av2 = [1./(zv2-zC), 1./(zv2-zC).^2, 1./(zv2-zrest.'), 1./(zv2-zout.'), (zv2/sc).^(0:pdeg)];
info.fit = max(abs(A*c-om));
info.val = max(abs(Av*c-conj(zv)));
info.val2 = max(abs(Av2*c-conj(zv2)));
end
