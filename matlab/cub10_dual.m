function cub10_dual()
%CUB10_DUAL   Error of the rule for f = 1/(z - xi), xi outside the ellipse.
%   With rho_in the interior-pole part of the fit and S the Schwarz function
%   of the ellipse, S(z) ~ kappa*z at infinity, kappa = (a-b)/(a+b),
%
%       I - I_n = pi*(rho_in(xi) - S(xi) + kappa*xi).
%
%   Checked at 18 points for degree 16. Also prints the conditioning of the
%   least-squares matrix. (Question of A. Horning.)
cfg = cub_config();
a = 1; b = 0.6; c = sqrt(a^2 - b^2);
t = 2*pi*(0:cfg.M-1)'/cfg.M;
zs = a*cos(t) + 1i*b*sin(t); om = conj(zs);

deg = 16;
w0 = warning('off', 'all');
[~, pol] = aaa(om, zs, 'degree', deg, 'sign', 1, 'lawson', 0, ...
               'tol', 0, 'cleanup', 0);
warning(w0);
pol = pol(isfinite(pol));
dmin = min(abs(zs - pol.'), [], 1);
pol = pol(dmin(:) > 1e-8*max(abs(zs - zs.'), [], 'all'));
in = inpolygon(real(pol), imag(pol), real(zs), imag(zs));
zk = pol(in); zout = pol(~in);
sc = max(abs(zs));
A = [1./(zs - zk.'), 1./(zs - zout.'), (zs/sc).^(0:cfg.pdeg)];
[co, lsinfo] = cub_lsq_svd(A, om);
ck = co(1:numel(zk)); wk = pi*ck;
epsf = max(abs(A*co - om));
nin = numel(zk);
Ain = A(:, 1:nin);
ninrm = sqrt(sum(abs(Ain).^2, 1));
sInterior = svd(Ain ./ ninrm, 'econ');
condInterior = sInterior(1) / sInterior(end);
if isempty(zout)
    maxExteriorPole = NaN;
    auxSpanResidual = NaN;
else
    [maxExteriorPole, jinf] = max(abs(zout));
    auxCol = A(:, nin + jinf);
    P = A(:, nin + numel(zout) + (1:cfg.pdeg+1));
    pNorm = sqrt(sum(abs(P).^2, 1));
    pAuxScaled = (P ./ pNorm) \ auxCol;
    pAux = pAuxScaled ./ pNorm.';
    auxSpanResidual = max(abs(P*pAux - auxCol)) / max(abs(auxCol));
end
fprintf('ellipse, deg %d: %d interior poles, fit residual = %.1e, ', ...
        deg, numel(zk), epsf);
condFull = lsinfo.sigma(1) / lsinfo.sigma(end);
fprintf('SVD rank = %d/%d, full cond. = %.1e, retained cond. = %.1e\n', ...
        lsinfo.rank, lsinfo.ncols, condFull, lsinfo.condEff);
fprintf('interior-pole cond. = %.1e, max exterior pole = %.1e, ', ...
        condInterior, maxExteriorPole);
fprintf('auxiliary-span residual = %.1e\n', auxSpanResidual);

% Schwarz function of the ellipse, cut on [-c,c], S(z) ~ kappa*z at infinity
S   = @(z) ((a^2 + b^2)*z - 2*a*b*z.*sqrt(1 - (c./z).^2))/c^2;
kap = (a - b)/(a + b);
zb  = a*cos(0.7) + 1i*b*sin(0.7);
assert(abs(S(zb) - conj(zb)) < 1e-14, 'Schwarz branch check failed');

zfun  = @(s) a*cos(s) + 1i*b*sin(s);
zpfun = @(s) -a*sin(s) + 1i*b*cos(s);
omfun = @(s) a*cos(s) - 1i*b*sin(s);

fprintf('%6s %6s %12s %12s %12s\n', ...
        'theta', 'scale', '|E(xi)|', '|pi*D(xi)|', '|E - pi*D|');
maxdef = 0; maxE = 0;
for u = [0.3, 1.2, 2.4]
    for s = [1.05 1.1 1.25 1.5 2 3]
        xi = s*(a*cos(u) + 1i*b*sin(u));
        f = @(z) 1./(z - xi);
        Iref = cub_ref(f, zfun, zpfun, omfun);
        In = sum(wk.*f(zk));
        E = Iref - In;
        D = sum(ck./(xi - zk)) - (S(xi) - kap*xi);
        maxdef = max(maxdef, abs(E - pi*D));
        maxE = max(maxE, abs(E));
        fprintf('%6.2f %6.2f %12.2e %12.2e %12.2e\n', ...
                u, s, abs(E), abs(pi*D), abs(E - pi*D));
    end
end
fprintf('identity defect: max |E - pi*D| = %.2e (max |E| = %.2e)\n', ...
        maxdef, maxE);
end

