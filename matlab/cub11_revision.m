function cub11_revision()
%CUB11_REVISION   Further measurements quoted in the paper.
%   (A) Ellipse and starfish: the same AAA poles fitted with four bases
%       (with or without the exterior poles, with or without the
%       polynomial). Error for f = exp(z) and Lambda_n.
%   (B) Ellipse: decay of the validation residual and of the error per
%       node against rho^(-2n), rho = (a+b)/c, and distance of the nodes
%       to the Gauss-Gegenbauer points.
%   (C) Square: the bound of Theorem 2.1 against the error, and Lambda_n.
%   (D) Ellipse: a harmonic integrand u = f + conj(g).
%   (E) Starfish: as (C).
cfg = cub_config();
M = cfg.M; L = cfg.pdeg;
t = 2*pi*(0:M-1)'/M;
f1 = @(z) exp(z);

%% (A) four bases, same poles
a = 1; b = 0.6; c = sqrt(a^2 - b^2);
zE = a*cos(t) + 1i*b*sin(t); omE = conj(zE);
[xj, wj] = jacpts(240, 0.5, 0.5);
IE = 2*a*b*(wj*f1(c*xj));
areaE = pi*a*b;
ablation('ellipse', zE, omE, IE, areaE, 4:4:44, L, f1);

rf  = @(s) 1 + 0.3*cos(5*s); rfp = @(s) -1.5*sin(5*s);
zf  = @(s) rf(s).*exp(1i*s); zfp = @(s) (rfp(s) + 1i*rf(s)).*exp(1i*s);
zS = zf(t); omS = conj(zS);
IS = cub_ref(f1, zf, zfp, @(s) conj(zf(s)));
areaS = real(cub_ref(@(z) 1 + 0*z, zf, zfp, @(s) conj(zf(s))));
ablation('starfish', zS, omS, IS, areaS, 30:30:150, L, f1);

%% (B) ellipse, rates
rho = (a + b)/c;
fprintf('\n(B) ellipse: rho = (a+b)/c = %.4f, rho^-2 = %.4f, log(rho^2) = %.4f per node\n', ...
        rho, rho^-2, log(rho^2));
tv = 2*pi*((0:4*M-1)' + 0.5)/(4*M);
zEv = a*cos(tv) + 1i*b*sin(tv); omEv = conj(zEv);
nn = []; ee = []; ev = []; dd = []; nbig = [];
fprintf('%4s %3s %4s %10s %10s %10s\n', 'deg', 'n', 'nbig', 'err', 'eps_val', 'dist2GG');
for deg = 2:2:44
    [zk, wk, epsv] = cub_rule(zE, omE, deg, L, zEv, omEv);
    if isempty(zk), continue; end
    big = abs(wk) > 1e-8*max(abs(wk));
    zb = zk(big);
    [xg, ~] = jacpts(numel(zb), 0.5, 0.5);
    dist = max(min(abs(zb - c*xg.'), [], 2));
    nn(end+1) = numel(zk); nbig(end+1) = numel(zb); %#ok<*AGROW>
    ee(end+1) = abs(sum(wk.*f1(zk)) - IE)/abs(IE); ev(end+1) = epsv; dd(end+1) = dist;
    fprintf('%4d %3d %4d %10.2e %10.2e %10.2e\n', deg, nn(end), nbig(end), ee(end), ev(end), dd(end));
end
msk = ev > 1e-13 & ev < 1e-2;
pe = polyfit(nbig(msk), log(ev(msk)), 1);
pr = polyfit(nbig(msk), log(ee(msk)), 1);
fprintf('fit over %d points: eps_val ~ exp(%.3f n_big), err ~ exp(%.3f n_big); rho^-2n = exp(%.3f n)\n', ...
        nnz(msk), pe(1), pr(1), -log(rho^2));
msk2 = ee > 1e-13 & ee < 1e-2;
pr2 = polyfit(nn(msk2), log(ee(msk2)), 1);
fprintf('err vs all interior poles n: exp(%.3f n)\n', pr2(1));

%% (C) square, bound against error
fprintf('\n(C) square [-1,1]^2\n');
mh = 240; dmin = 1e-12;
lg = linspace(log(dmin), 0, mh)';
zq = square_samples(exp(lg)); h = lg(2) - lg(1);
zqv = square_samples(exp(lg(1:end-1) + h/2));
omq = conj(zq); omqv = conj(zqv);
I1 = (exp(1) - exp(-1))*2*sin(1);
glen = 8; fmax = exp(1);
fprintf('%4s %4s %10s %10s %10s %10s %8s\n', 'deg', 'n', 'err', 'eps_val', 'bound', 'ratio', 'Lambda');
for deg = 20:20:220
    [zk, wk, epsv] = cub_rule(zq, omq, deg, L, zqv, omqv);
    if isempty(zk), continue; end
    err = abs(sum(wk.*f1(zk)) - I1);
    bnd = 0.5*glen*epsv*fmax;
    fprintf('%4d %4d %10.2e %10.2e %10.2e %10.2e %8.3f   |sum w - 4| = %.1e\n', ...
            deg, numel(zk), err/abs(I1), epsv, bnd/abs(I1), err/bnd, sum(abs(wk))/4, abs(sum(wk) - 4));
end

%% (D) harmonic integrand
fprintf('\n(D) harmonic integrand u = f + conj(g), ellipse, f = exp(z), g = cos(z)\n');
g = @(z) cos(z);
IEg = 2*a*b*(wj*g(c*xj));
for deg = [20 30 44]
    [zk, wk] = cub_rule(zE, omE, deg, L);
    Iu = sum(wk.*f1(zk)) + conj(sum(wk.*g(zk)));
    fprintf('  deg %2d: |I_n(u) - I(u)| / |I(u)| = %.2e\n', deg, abs(Iu - (IE + conj(IEg)))/abs(IE + conj(IEg)));
end

%% (E) starfish, bound against error
fprintf('\n(E) starfish: bound vs error\n');
tvS = 2*pi*((0:4*M-1)' + 0.5)/(4*M);
zSv = zf(tvS); omSv = conj(zSv);
glenS = (2*pi/numel(tvS))*sum(abs(zfp(tvS)));
fmaxS = max(abs(f1(zS)));
fprintf('%4s %3s %10s %10s %10s %10s %8s\n', 'deg', 'n', 'err', 'eps_val', 'bound', 'ratio', 'Lambda');
for deg = 30:30:150
    [zk, wk, epsv] = cub_rule(zS, omS, deg, L, zSv, omSv);
    err = abs(sum(wk.*f1(zk)) - IS);
    bnd = 0.5*glenS*epsv*fmaxS;
    fprintf('%4d %3d %10.2e %10.2e %10.2e %10.2e %8.3f\n', deg, numel(zk), ...
            err/abs(IS), epsv, bnd/abs(IS), err/bnd, sum(abs(wk))/areaS);
end
end

function ablation(name, zs, om, Iref, area, degs, L, f1)
fprintf('\n(A) %s: error / Lambda_n for four least-squares bases (same AAA poles)\n', name);
fprintf('%4s %3s | %9s %6s | %9s %6s | %9s %6s | %9s %6s\n', 'deg', 'n', ...
        'full', 'Lam', 'int only', 'Lam', 'int+poly', 'Lam', 'int+ext', 'Lam');
for deg = degs
    [zk, zout] = poles(zs, om, deg);
    if isempty(zk), continue; end
    out = zeros(1, 8);
    opts = {[true true], [false false], [false true], [true false]};
    for v = 1:4
        [wk, ~] = fit(zs, om, zk, zout, L, opts{v}(1), opts{v}(2));
        out(2*v-1) = abs(sum(wk.*f1(zk)) - Iref)/abs(Iref);
        out(2*v) = sum(abs(wk))/area;
    end
    fprintf('%4d %3d | %9.2e %6.3f | %9.2e %6.3f | %9.2e %6.3f | %9.2e %6.3f\n', ...
            deg, numel(zk), out);
end
end

function [zk, zout] = poles(zs, om, deg)
ws = warning('off', 'all');
[~, pol] = aaa(om, zs, 'degree', deg, 'sign', 1, 'lawson', 0, 'tol', 0, 'cleanup', 0);
warning(ws);
pol = pol(isfinite(pol));
dmin = min(abs(zs - pol.'), [], 1);
pol = pol(dmin(:) > 1e-8*max(abs(zs - zs.'), [], 'all'));
in = inpolygon(real(pol), imag(pol), real(zs), imag(zs));
zk = pol(in); zout = pol(~in);
end

function [wk, epsf] = fit(zs, om, zk, zout, L, useExt, usePoly)
sc = max(abs(zs));
A = 1./(zs - zk.');
if useExt, A = [A, 1./(zs - zout.')]; end
if usePoly, A = [A, (zs/sc).^(0:L)]; end
c = cub_lsq_svd(A, om);
epsf = max(abs(A*c - om));
wk = pi*c(1:numel(zk));
end

function zs = square_samples(d)
d = d(:);
u = unique([-1 + d; 1 - d]);
zs = [u - 1i; 1 + 1i*u; -u + 1i; -1 - 1i*u];
end
