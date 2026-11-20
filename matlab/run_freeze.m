 function run_freeze()
%RUN_FREEZE   Run cub1 to cub10 in order and print the time each one takes.
fprintf('=== FREEZE RUN %s ===\n', datestr(now));
fprintf('MATLAB %s\n', version);
fns = {@cub1_disk, @cub2_ellipse, @cub3_domains, @cub3b_qd, ...
       @cub3c_resolution, @cub4_potential, @cub4b_resolution, ...
       @cub5_kernels, @cub6_perturbed, @cub7_details, ...
       @cub8_square, @cub9_neumann, @cub10_dual};
for k = 1:numel(fns)
    nm = func2str(fns{k});
    fprintf('\n========== %s ==========\n', nm);
    t0 = tic;
    fns{k}();
    fprintf('---- %s done in %.1f s\n', nm, toc(t0));
    close all
end
fprintf('=== FREEZE COMPLETE ===\n');
end
