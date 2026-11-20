function cfg = cub_config()
%CUB_CONFIG   Parameters and plot defaults shared by the cub* scripts.
%   Chebfun (jacpts, legpts, aaa) has to be on the path. If it is not, put
%   the addpath calls in local_paths.m, see local_paths_template.m, or set
%   the environment variable CHEBFUN_PATH.
here = fileparts(mfilename('fullpath'));
cfg.here = here;
if exist(fullfile(here, 'local_paths.m'), 'file') == 2
    run(fullfile(here, 'local_paths.m'));
end
if exist('aaa', 'file') ~= 2 && ~isempty(getenv('CHEBFUN_PATH'))
    addpath(getenv('CHEBFUN_PATH'));
end
if exist('aaa', 'file') ~= 2 || exist('jacpts', 'file') ~= 2
    error('cub_config:chebfun', ...
          'Chebfun not found, see local_paths_template.m.');
end
cfg.deps.aaa = which('aaa');
cfg.deps.jacpts = which('jacpts');

cfg.M = 800;            % boundary points
cfg.pdeg = 10;          % polynomial degree in the least-squares fit

% plots
cfg.FS = 18;
cfg.col.node  = [0.85 0.10 0.10];
cfg.col.curve = [0.00 0.35 0.72];
cfg.col.gauss = [0.15 0.15 0.15];
cfg.col.alt   = [0.00 0.35 0.72];
cfg.col.violet = [0.45 0.20 0.55];
cfg.fig.position = [100 100 1000 380];
cfg.plot.curveWidth = 1.1;
cfg.plot.markerSize = 6.5;
cfg.plot.nodeSize = 12;
cfg.plot.errorLimits = [1e-16 10];
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultFigureRenderer', 'painters');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultAxesLineWidth', 0.6);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesGridAlpha', 0.08);
set(groot, 'defaultAxesMinorGridAlpha', 0.04);
set(groot, 'defaultAxesBox', 'on');
end
