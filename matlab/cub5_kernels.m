function cub5_kernels()
%CUB5_KERNELS   Cauchy kernel in place of the logarithm.
%
%       I(z0) = int int_Omega f(z)/(z - z0) dA,    f = exp(z).
%
%   Boundary data W = (zbar - conj(z0))/(z - z0), which is bounded at z0.
%   Starfish with M = 1600 points, one target inside and one outside.
cfg = cub_config();
rf  = @(s) 1 + 0.3*cos(5*s);
rfp = @(s) -1.5*sin(5*s);
zf  = @(s) rf(s).*exp(1i*s);
zfp = @(s) (rfp(s) + 1i*rf(s)).*exp(1i*s);
f1 = @(z) exp(z);

tdir = 0.55;
zb = zf(tdir); nrm = 1i*zfp(tdir)/abs(zfp(tdir));
if real(conj(zb)*nrm) < 0, nrm = -nrm; end
cases = {0.20 + 0.10i, 'interior, central'; ...
         zb + 0.05*nrm, 'exterior, delta = 0.05'};

M = 1600;
t = 2*pi*(0:M-1)'/M;
zs = zf(t);
for j = 1:size(cases, 1)
    z0 = cases{j,1};
    Wf = @(s) (conj(zf(s)) - conj(z0))./(zf(s) - z0);
    om = Wf(t);
    Iref = cub_ref(f1, zf, zfp, Wf);
    best = inf; bn = 0;
    for deg = 10:10:120
        [zk, wk, ~] = cub_rule(zs, om, deg, cfg.pdeg);
        if isempty(zk), continue; end
        e = abs(sum(wk.*f1(zk)) - Iref)/abs(Iref);
        if e < best, best = e; bn = numel(zk); end
    end
    fprintf('Cauchy kernel, %-22s I = %.12f%+.12fi, best rel err %.2e at n = %d\n', ...
            cases{j,2}, real(Iref), imag(Iref), best, bn);
end
end
