function run_all_cub()
%RUN_ALL_CUB   Short run: cub1 to cub4.
fns = {@cub1_disk, @cub2_ellipse, @cub3_domains, @cub4_potential};
for k = 1:numel(fns)
    fprintf('\n%s\n', func2str(fns{k}));
    tic, fns{k}(); toc
    close all
end
end
