function I = cub_ref(fh, zfun, zpfun, omfun)
%CUB_REF   Reference value of an area integral.
%   I = cub_ref(f, z, zp, W) is (1/2i) * int f(z(t)) W(t) z'(t) dt over
%   [0, 2*pi] by the trapezoidal rule with 2^14 points, checked against
%   2^13 points. f, z, zp, W are function handles; z, zp, W take t.
I1 = trap(2^13); I2 = trap(2^14);
assert(abs(I1 - I2) < 1e-12 * max(1, abs(I2)), ...
       'reference self-check failed: %.2e', abs(I1 - I2));
I = I2;
    function I = trap(M)
        t = 2*pi*(0:M-1)'/M;
        I = (1/(2i)) * (2*pi/M) * sum(fh(zfun(t)) .* omfun(t) .* zpfun(t));
    end
end
