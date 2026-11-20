function cub9_neumann()
%CUB9_NEUMANN   Neumann oval (Figure 5).
%   Omega = phi(D), phi(zeta) = zeta/(1 - q^2*zeta^2), q = 1/2, is a
%   quadrature domain with two simple nodes:
%
%       int int f dA = pi*c*(f(x0) + f(-x0)),
%       x0 = q/(1 - q^4) = 8/15,    c = (1 + q^4)/(2*(1 - q^4)^2) = 136/225.
%
%   The script compares the computed nodes and weights with x0 and pi*c.
%   In cub3b_qd the two nodes coincide.
cfg = cub_config();
FS = cfg.FS;
q = 0.5;
phi  = @(zeta) zeta./(1 - q^2*zeta.^2);
phip = @(zeta) (1 + q^2*zeta.^2)./(1 - q^2*zeta.^2).^2;
t = 2*pi*(0:cfg.M-1)'/cfg.M;
zs = phi(exp(1i*t)); om = conj(zs);
xs  = q/(1 - q^4);
cex = (1 + q^4)/(2*(1 - q^4)^2);
fprintf('exact identity: nodes +-%.15f, weight pi * %.15f\n', xs, cex);

% exact rule against the boundary integral
f1 = @(z) exp(z); f2 = @(z) 1./(1.6 - z);
fs = {f1, f2, @(z) z.^2, @(z) z.^7, @(z) cos(3*z)};
zfun  = @(s) phi(exp(1i*s));
zpfun = @(s) 1i*exp(1i*s).*phip(exp(1i*s));
omfun = @(s) conj(phi(exp(1i*s)));
Iex = @(f) pi*cex*(f(xs) + f(-xs));
for k = 1:numel(fs)
    Ik = cub_ref(fs{k}, zfun, zpfun, omfun);
    fprintf('  identity vs boundary ref, f%d: %.2e\n', k, ...
            abs(Ik - Iex(fs{k}))/max(1, abs(Ik)));
end

% AAA degrees 4:2:24, keep the most accurate rule
dg = []; an = []; ea = []; best = {[], [], inf, 0};
for deg = 4:2:24
    [zkd, wkd, epsd] = cub_rule(zs, om, deg, cfg.pdeg);
    if isempty(zkd), continue; end
    e = 0;
    for k = 1:numel(fs)
        Ik = Iex(fs{k});
        e = max(e, abs(sum(wkd.*fs{k}(zkd)) - Ik)/max(1, abs(Ik)));
    end
    dg(end+1) = deg; an(end+1) = numel(zkd); ea(end+1) = e; %#ok<*AGROW>
    if e < best{3}, best = {zkd, wkd, e, epsd}; end
end
fprintf('interior node count across degrees %d..%d: min %d, max %d\n', ...
        dg(1), dg(end), min(an), max(an));
[zk, wk, epsf] = deal(best{1}, best{2}, best{4});
[~, ord] = sort(abs(wk), 'descend'); zk = zk(ord); wk = wk(ord);
big = abs(wk) > 1e-10;
fprintf('computed rule: %d nodes carry weight (eps = %.1e)\n', sum(big), epsf);
for k = find(big).'
    fprintf('  node %+.15f %+.2ei   weight/pi = %.15f %+.2ei\n', ...
            real(zk(k)), imag(zk(k)), real(wk(k))/pi, imag(wk(k))/pi);
end
zkb = zk(big); wkb = wk(big);
[~, ip] = max(real(zkb)); [~, im] = min(real(zkb));
fprintf('  node errors:   %.2e  %.2e\n', abs(zkb(ip) - xs), abs(zkb(im) + xs));
fprintf('  weight errors: %.2e  %.2e\n', ...
        abs(wkb(ip) - pi*cex), abs(wkb(im) - pi*cex));
fprintf('  node separation: %.6f (exact %.6f)\n', ...
        abs(zkb(ip) - zkb(im)), 2*xs);
maxe = 0;
for k = 1:numel(fs)
    Ik = Iex(fs{k});
    maxe = max(maxe, abs(sum(wk.*fs{k}(zk)) - Ik)/max(1, abs(Ik)));
end
fprintf('  exactness over 5 analytic integrands: max scaled err = %.2e\n', maxe);
fprintf('  sum of weights - exact area 2*pi*c: %.2e\n', abs(sum(wk) - 2*pi*cex));

% figure
figure('Color','w','Position',cfg.fig.position);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile
plot(real(zs([1:end 1])), imag(zs([1:end 1])), '-', 'Color', cfg.col.curve, ...
     'LineWidth', cfg.plot.curveWidth); hold on
plot(real(zkb), imag(zkb), '.', 'Color', cfg.col.node, 'MarkerSize', 14)
plot([xs -xs], [0 0], 'k+', 'MarkerSize', 9, 'LineWidth', 1.1)
axis equal; box on
xlabel('Re $z$', 'FontSize', FS); ylabel('Im $z$', 'FontSize', FS)
nexttile
semilogy(dg, ea, 'o', 'Color', cfg.col.node, ...
         'MarkerSize', cfg.plot.markerSize, 'LineWidth', 1)
ylim(cfg.plot.errorLimits); grid on; box on; set(gca,'YMinorGrid','off')
xlabel('AAA degree', 'FontSize', FS)
ylabel('max scaled error', 'FontSize', FS)
for ax = findall(gcf, 'Type', 'axes')'
    ax.Toolbar.Visible = 'off';
end
exportgraphics(gcf, fullfile(cfg.here, 'cub9_neumann.pdf'), ...
               'ContentType', 'vector', 'BackgroundColor', 'white');
end
