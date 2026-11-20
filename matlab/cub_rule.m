function [zk, wk, epsf, zout] = cub_rule(zs, om, deg, pdeg, zsv, omv)
%CUB_RULE   Cubature rule from boundary data.
%   [zk, wk, epsf, zout] = cub_rule(zs, om, deg, pdeg)
%
%   zs    points on the boundary curve, in order (column)
%   om    boundary data W at zs, where dW/dzbar = w; conj(zs) if w = 1
%   deg   AAA degree
%   pdeg  degree of the polynomial part of the fit
%
%   zk    nodes, the AAA poles inside the curve
%   wk    weights, pi times the fitted residues
%   epsf  max |rho - om| on zs, or on zsv for cub_rule(..., zsv, omv)
%   zout  poles outside the curve, kept in the fit but without weight
%
%   sum(wk.*f(zk)) approximates the integral of f*w over the domain, with
%   error at most (length/2) * max|rho - W| * max|f| on the curve.
if nargin < 5, zsv = []; omv = []; end

% poles from AAA, no Lawson steps, no cleanup
w0 = warning('off', 'all');
[~, pol] = aaa(om, zs, 'degree', deg, 'sign', 1, 'lawson', 0, ...
               'tol', 0, 'cleanup', 0);
warning(w0);

% drop poles at infinity or within 1e-8 diameters of a sample point
pol = pol(isfinite(pol));
diam = max(abs(zs - zs.'), [], 'all');
dmin = min(abs(zs - pol.'), [], 1);
pol = pol(dmin(:) > 1e-8*diam);

in = inpolygon(real(pol), imag(pol), real(zs), imag(zs));
zk = pol(in);
zout = pol(~in);

% least squares on [1/(z-zk), 1/(z-zout), (z/sc)^j]
sc = max(abs(zs));
A = [1./(zs - zk.'), 1./(zs - zout.'), (zs/sc).^(0:pdeg)];
c = cub_lsq_svd(A, om);
if isempty(zsv)
    epsf = max(abs(A*c - om));
else
    Av = [1./(zsv - zk.'), 1./(zsv - zout.'), (zsv/sc).^(0:pdeg)];
    epsf = max(abs(Av*c - omv));
end
wk = pi*c(1:numel(zk));
end
