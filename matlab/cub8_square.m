function cub8_square()
%CUB8_SQUARE   The square [-1,1]^2 (Figure 6).
%   zbar is linear in z on each edge, with a different linear function on
%   each, so the boundary data has a kink at every corner. The boundary
%   points are clustered exponentially at the corners, down to 1e-12.
%   Prints the errors for f = exp(z) and f = 1/(1.6 - z), the distance of
%   the nodes to the diagonals and to the corners, and a fit of the error
%   to exp(-C*sqrt(n)).
cfg = cub_config();
FS = cfg.FS;

% boundary points, clustered at both ends of each edge
mh   = 240;                        % samples per half-edge
dmin = 1e-12;                      % smallest distance to a corner
lg  = linspace(log(dmin), 0, mh)';
zs  = square_samples(exp(lg));
h   = lg(2) - lg(1);
zsv = square_samples(exp(lg(1:end-1) + h/2));   % offset validation grid
om  = conj(zs); omv = conj(zsv);
fprintf('square: %d boundary samples, %d validation samples\n', ...
        numel(zs), numel(zsv));

% references: closed form for exp, tensor Gauss-Legendre for both
f1 = @(z) exp(z);
f2 = @(z) 1./(1.6 - z);
I1 = (exp(1) - exp(-1)) * 2*sin(1);             % int exp(x) dx * int exp(iy) dy
[xg, wg] = legpts(120);
[X, Y] = ndgrid(xg, xg);
I1g = wg * f1(X + 1i*Y) * wg';
I2  = wg * f2(X + 1i*Y) * wg';
fprintf('square: tensor GL vs closed form for exp: %.2e\n', abs(I1g - I1));

an = []; e1 = []; e2 = []; ee = [];
keepin = []; keepout = []; keepw = [];
for deg = 20:20:220
    [zk, wk, epsf, zout] = cub_rule(zs, om, deg, cfg.pdeg, zsv, omv);
    if isempty(zk), continue; end
    an(end+1) = numel(zk); %#ok<*AGROW>
    e1(end+1) = abs(sum(wk.*f1(zk)) - I1)/abs(I1);
    e2(end+1) = abs(sum(wk.*f2(zk)) - I2)/abs(I2);
    ee(end+1) = epsf;
    fprintf(['deg %3d: n = %3d, err(exp) = %.2e, err(pole) = %.2e, ' ...
             'eps = %.1e, |sum w - 4| = %.2e\n'], ...
            deg, numel(zk), e1(end), e2(end), epsf, abs(sum(wk) - 4));
    if isempty(keepin) || e1(end) <= min(e1)
        keepin = zk; keepout = zout; keepw = wk;
    end
end
[be, bi] = min(e1);
fprintf('square: best rel err %.2e at n = %d\n', be, an(bi));

% distance of the nodes to the diagonals and to the nearest corner
wbig = abs(keepw) > 1e-8 * max(abs(keepw));
zkb  = keepin(wbig);
ddiag = min(abs(real(zkb) - imag(zkb)), abs(real(zkb) + imag(zkb)))/sqrt(2);
corners = [1+1i; -1+1i; -1-1i; 1-1i];
dcorn = min(abs(zkb - corners.'), [], 2);
fprintf(['nodes carrying weight: %d;  dist to diagonals: median = %.2e, ' ...
         'max = %.2e;  min dist to a corner = %.2e\n'], ...
        numel(zkb), median(ddiag), max(ddiag), min(dcorn));

% fit log(err) against sqrt(n) where 1e-13 < err < 1e-1
msk = e1 < 1e-1 & e1 > 1e-13;
p = polyfit(sqrt(an(msk)), log(e1(msk)), 1);
fprintf('root-exponential fit: err ~ exp(%+.2f sqrt(n))\n', p(1));

% figure
figure('Color','w','Position',cfg.fig.position);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile
plot([1 -1 -1 1 1], [1 1 -1 -1 1], '-', 'Color', cfg.col.curve, ...
     'LineWidth', cfg.plot.curveWidth); hold on
plot([-1 1], [-1 1], ':', 'Color', [.6 .6 .6], 'LineWidth', 0.8)
plot([-1 1], [1 -1], ':', 'Color', [.6 .6 .6], 'LineWidth', 0.8)
plot(real(keepout), imag(keepout), '.', 'Color', [.65 .65 .65], 'MarkerSize', 8)
plot(real(zkb), imag(zkb), '.', 'Color', cfg.col.node, ...
     'MarkerSize', cfg.plot.nodeSize)
axis equal; box on; xlim([-2.1 2.1]); ylim([-2.1 2.1])
xlabel('Re $z$', 'FontSize', FS); ylabel('Im $z$', 'FontSize', FS)
nexttile
sq = sqrt(an);
semilogy(sq, e1, 'o', 'Color', cfg.col.node, ...
         'MarkerSize', cfg.plot.markerSize, 'LineWidth', 1); hold on
semilogy(sq, e2, 's', 'Color', cfg.col.alt, ...
         'MarkerSize', cfg.plot.markerSize)
sfit = linspace(min(sq), max(sq), 100);
semilogy(sfit, exp(p(2) + p(1)*sfit), '--', 'Color', [.45 .45 .45], 'LineWidth', 1)
ylim(cfg.plot.errorLimits); grid on; box on; set(gca,'YMinorGrid','off')
xlabel('$\sqrt{n}$', 'FontSize', FS)
ylabel('relative error', 'FontSize', FS)
text(sq(2), 3e-2, '$f = e^z$', 'Color', cfg.col.node, 'FontSize', FS-1)
text(sq(4), 1e-9, '$f = 1/(1.6-z)$', 'Color', cfg.col.alt, 'FontSize', FS-1)
for ax = findall(gcf, 'Type', 'axes')'
    ax.Toolbar.Visible = 'off';
end
exportgraphics(gcf, fullfile(cfg.here, 'cub8_square.pdf'), ...
               'ContentType', 'vector', 'BackgroundColor', 'white');
end

function zs = square_samples(d)
% Points on the boundary of [-1,1]^2 at distances d (in (0,1]) from each
% corner, counterclockwise starting near -1-1i.
d = d(:);
u = unique([-1 + d; 1 - d]);          % coordinate along an edge
zs = [u - 1i;                         % bottom:  -1-1i -> 1-1i
      1 + 1i*u;                       % right:    1-1i -> 1+1i
      -u + 1i;                        % top:      1+1i -> -1+1i
      -1 - 1i*u];                     % left:    -1+1i -> -1-1i
end
