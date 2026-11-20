function run_reproduce()
%RUN_REPRODUCE   All computations of the paper.
%   Calls run_freeze, then cub11_revision and both modes of cub12_checks,
%   and copies the seven figure PDFs to ../figures.
here = fileparts(mfilename('fullpath'));
run_freeze();
cub11_revision();
cub12_checks('basis');
cub12_checks('ellipse');
close all

% Figures 1 to 7
figs = {'cub2_ellipse', 'cub6_perturbed', 'cub3_domains', 'cub3b_qd', ...
        'cub9_neumann', 'cub8_square', 'cub4_potential'};
figdir = fullfile(here, '..', 'figures');
if ~exist(figdir, 'dir'), mkdir(figdir); end
for k = 1:numel(figs)
    copyfile(fullfile(here, [figs{k} '.pdf']), figdir);
end
end
