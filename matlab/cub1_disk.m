function cub1_disk()
%CUB1_DISK   Disk: the rule reduces to the mean value property.
%   On |z - z0| = R we have zbar = conj(z0) + R^2/(z - z0), so the boundary
%   data have one pole, at the centre, with residue R^2. The script checks
%   that the node of largest weight is z0 with weight pi*R^2, for the unit
%   disk and for the disk of radius 0.7 about 0.5 + 0.2i.
cfg = cub_config();
t = 2*pi*(0:cfg.M-1)'/cfg.M;

fprintf('--- unit disk ---\n');
zs = exp(1i*t);
[zk, wk, epsf] = cub_rule(zs, conj(zs), 12, cfg.pdeg);
[~, km] = max(abs(wk));
fprintf('nodes: %d;  dominant node = %.2e%+.2ei;  weight - pi = %.2e;  eps = %.1e\n', ...
        numel(zk), real(zk(km)), imag(zk(km)), abs(wk(km) - pi), epsf);
fs = {@(z) exp(z), @(z) 1./(3 - z), @(z) cos(2*z)};
for k = 1:3
    In = sum(wk .* fs{k}(zk));
    Iref = pi * fs{k}(0);                       % mean value property
    fprintf('  f%d: error = %.2e\n', k, abs(In - Iref));
end

fprintf('--- disk of radius 0.7 centered at 0.5 + 0.2i ---\n');
z0 = 0.5 + 0.2i; R = 0.7;
zs = z0 + R*exp(1i*t);
[zk, wk, epsf] = cub_rule(zs, conj(zs), 12, cfg.pdeg);
[~, kmain] = max(abs(wk));
fprintf('nodes: %d;  main node - center = %.2e;  weight - pi R^2 = %.2e;  eps = %.1e\n', ...
        numel(zk), abs(zk(kmain) - z0), abs(wk(kmain) - pi*R^2), epsf);
In = sum(wk .* exp(zk));
fprintf('  exp: error vs pi R^2 exp(center) = %.2e\n', abs(In - pi*R^2*exp(z0)));
end
