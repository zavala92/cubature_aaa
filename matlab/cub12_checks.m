function cub12_checks(mode)
%CUB12_CHECKS   Table 1 and the rate on the ellipse.
%   cub12_checks('basis')    the four fits of Table 1: ellipse, degrees 24
%                            and 44; starfish, degree 150
%   cub12_checks('ellipse')  decay of the validation residual on the
%                            ellipse, and the Gauss-Gegenbauer remainder
%                            for n = 3, 7, 11
cfg=cub_config(); t=2*pi*(0:cfg.M-1)'/cfg.M;
fprintf('AAA implementation: %s\n',cfg.deps.aaa);
if strcmp(mode,'basis')
    a=1;b=.6;zs=a*cos(t)+1i*b*sin(t);
    [x,w]=jacpts(240,.5,.5); I=2*a*b*(w*exp(sqrt(a*a-b*b)*x));
    basis('ellipse',zs,I,pi*a*b,[24,44],cfg.pdeg);
    zf=@(s) (1+.3*cos(5*s)).*exp(1i*s);
    zfp=@(s) (-1.5*sin(5*s)+1i*(1+.3*cos(5*s))).*exp(1i*s);
    I=cub_ref(@exp,zf,zfp,@(s) conj(zf(s)));
    basis('starfish',zf(t),I,pi*(1+.3^2/2),150,cfg.pdeg);
elseif strcmp(mode,'ellipse')
    a=1;b=.6;c=sqrt(a*a-b*b); zs=a*cos(t)+1i*b*sin(t);
    tv=2*pi*((0:4*cfg.M-1)'+.5)/(4*cfg.M);
    zv=a*cos(tv)+1i*b*sin(tv); kap=(a-b)/(a+b);
    ww=(zv+zv.*sqrt(1-(c./zv).^2))/c;
    H=2*a*b./(c*ww); res=[]; nodes=[]; degrees=[];
    for deg=2:2:26
        [zk,~,ev]=cub_rule(zs,conj(zs),deg,cfg.pdeg,zv,conj(zv));
        res(end+1)=ev;nodes(end+1)=numel(zk);degrees(end+1)=deg;
    end
    use=res>1e-13 & res<1e-2;
    pf=polyfit(nodes(use),log(res(use)),1);
    fprintf('AAA slope %.6f; fitted degrees %s; nodes %s\n',pf(1),mat2str(degrees(use)),mat2str(nodes(use)));
    for n=[3,7,11]
        [x,w]=jacpts(n,.5,.5);
        Hn=(2*a*b/pi)*sum(w(:).'./(zv-c*x(:).'),2);
        errTheory=2*a*b/c*ww.^(-2*n-1).*(1-ww.^(-2))./(1-ww.^(-2*n-2));
        fprintf('GG n%d: boundary residual %.3e; explicit remainder defect %.3e\n',n, ...
            max(abs(kap*zv+Hn-conj(zv))),max(abs((H-Hn)-errTheory)));
        assert(max(abs((H-Hn)-errTheory))<1e-12);
    end
else
    error('mode must be basis or ellipse');
end
end
function basis(name,zs,I,area,degrees,L)
for deg=degrees
    [zk,~,~,zout]=cub_rule(zs,conj(zs),deg,L);
    sc=max(abs(zs));
    opts=[0 0;0 1;1 0;1 1];
    labels={'interior','interior+poly','interior+exterior','full'};
    for j=1:4
        A=1./(zs-zk.');
        if opts(j,1),A=[A,1./(zs-zout.')];end
        if opts(j,2),A=[A,(zs/sc).^(0:L)];end
        [co,info]=cub_lsq_svd(A,conj(zs));
        As=A./sqrt(sum(abs(A).^2,1));
        assert(info.rank==rank(As,info.tol));
        assert(info.tol==max(size(A))*eps(info.sigma(1)));
        wk=pi*co(1:numel(zk));
        fprintf('%s degree %d n%d %s: cubature %.6e Lambda %.9g fit_residual %.6e rank %d/%d\n', ...
            name,deg,numel(zk),labels{j},abs(sum(wk.*exp(zk))-I)/abs(I), ...
            sum(abs(wk))/area,max(abs(A*co-conj(zs))),info.rank,size(A,2));
    end
end
end
