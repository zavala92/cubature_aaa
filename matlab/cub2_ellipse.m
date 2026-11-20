function cub2_ellipse()
%CUB2_ELLIPSE   Ellipse with semi-axes a = 1, b = 0.6 (Figure 1).
%   The area integral reduces to the focal segment [-c, c], c^2 = a^2 - b^2:
%
%       I = (2ab/c^2) * int_{-c}^{c} f(x) sqrt(c^2 - x^2) dx,
%
%   so Gauss-Gegenbauer quadrature on [-c, c] is the rule to compare with.
%   Left panel: nodes of the first rule with n >= 14 and the Gauss points.
%   Right panel: errors for f = exp(z) against the number of nodes.
cfg = cub_config();
a = 1; b = 0.6; c = sqrt(a^2 - b^2);
t = 2*pi*(0:cfg.M-1)'/cfg.M;
zs = a*cos(t) + 1i*b*sin(t);
om = conj(zs);

f1 = @(z) exp(z);
f2 = @(z) 1./(1.5 - z);
[xj, wj] = jacpts(240, 0.5, 0.5);           % reference on the focal segment
I1 = 2*a*b * (wj * f1(c*xj));
I2 = 2*a*b * (wj * f2(c*xj));
Ib = cub_ref(f1, @(s) a*cos(s)+1i*b*sin(s), @(s) -a*sin(s)+1i*b*cos(s), ...
             @(s) a*cos(s)-1i*b*sin(s));
assert(abs(I1 - Ib) < 1e-12 * abs(I1), 'focal-segment identity check failed');
fprintf('ellipse: focal identity confirmed to %.1e; I(exp) = %.15f\n', ...
        abs(I1 - Ib), real(I1));

an = []; e1 = []; e2 = []; keep = [];
for deg = 2:2:44
    [zk, wk, ~] = cub_rule(zs, om, deg, cfg.pdeg);
    if isempty(zk), continue; end
    an(end+1) = numel(zk); %#ok<*AGROW>
    e1(end+1) = abs(sum(wk.*f1(zk)) - I1)/abs(I1);
    e2(end+1) = abs(sum(wk.*f2(zk)) - I2)/abs(I2);
    if numel(zk) >= 14 && isempty(keep), keep = zk; end
end
% Gauss-Gegenbauer with 1 to 20 points
gn = 1:20; g1 = zeros(size(gn)); g2 = g1;
for k = gn
    [xg, wg] = jacpts(k, 0.5, 0.5);
    g1(k) = abs(2*a*b*(wg*f1(c*xg)) - I1)/abs(I1);
    g2(k) = abs(2*a*b*(wg*f2(c*xg)) - I2)/abs(I2);
end
fprintf('best AAA rule: n = %d, rel err = %.2e (exp)\n', an(end), e1(end));
fprintf('max |Im node| (should be ~0, nodes on focal segment): %.2e\n', ...
        max(abs(imag(keep))));

FS = cfg.FS;
figure('Color','w','Position',cfg.fig.position);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile
plot(real(zs([1:end 1])), imag(zs([1:end 1])), '-', 'Color', cfg.col.curve, ...
     'LineWidth', cfg.plot.curveWidth); hold on
plot([-c c], [0 0], 'k-', 'LineWidth', 0.8)
[xg, ~] = jacpts(numel(keep), 0.5, 0.5);
plot(c*xg, 0*xg, 'o', 'Color', cfg.col.gauss, 'MarkerSize', 7)
plot(real(keep), imag(keep), '.', 'Color', cfg.col.node, 'MarkerSize', 12)
plot([-c c], [0 0], 'k+', 'MarkerSize', 8, 'LineWidth', 1.1)
axis equal; box on
xlabel('Re $z$', 'FontSize', FS); ylabel('Im $z$', 'FontSize', FS)
nexttile
semilogy(an, e1, 'o', 'Color', cfg.col.node, ...
         'MarkerSize', cfg.plot.markerSize, 'LineWidth', 1); hold on
semilogy(gn, g1, 's', 'Color', cfg.col.gauss, ...
         'MarkerSize', cfg.plot.markerSize)
ylim(cfg.plot.errorLimits); grid on; box on; set(gca,'YMinorGrid','off')
xlabel('number of nodes $n$', 'FontSize', FS)
ylabel('relative error', 'FontSize', FS)
legend({'AAA cubature', 'Gauss--Gegenbauer'}, 'FontSize', FS-3, ...
       'Box', 'off', 'Location', 'northeast')
for ax = findall(gcf, 'Type', 'axes')'
    ax.Toolbar.Visible = 'off';
end
exportgraphics(gcf, fullfile(cfg.here, 'cub2_ellipse.pdf'), ...
               'ContentType', 'vector', 'BackgroundColor', 'white');
end
