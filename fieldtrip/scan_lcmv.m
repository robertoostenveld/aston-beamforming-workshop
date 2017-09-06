% LCMV scan

%%

cfg  = [];
cfg.headmodel = headmodel;
cfg.grid = sourcemodel;
cfg.method = 'lcmv';
cfg.channel = 'megplanar';
%cfg.latency = 0.060;
source_lcmv_active   = ft_sourceanalysis(cfg, timelock_active);

%cfg.latency = 0.060+0.150;
source_lcmv_baseline = ft_sourceanalysis(cfg, timelock_baseline);

%%

source_lcmv_active.time   = 0;
source_lcmv_baseline.time = 0;

cfg = [];
cfg.parameter = 'pow';
cfg.operation = 'log10(x1/x2)';
source_lcmv_relative = ft_math(cfg, source_lcmv_active, source_lcmv_baseline);

%%

cfg = [];
cfg.funparameter = 'pow';
cfg.location = 'max';
ft_sourceplot(cfg, source_lcmv_relative);

%%

cfg = [];
cfg.downsample = 2; % for memory reasons
cfg.parameter = 'pow';
source_lcmv_relative_int = ft_sourceinterpolate(cfg, source_lcmv_relative, mri);

%%
% PLEASE NOTE THAT i IS COMING FROM THE CALLING SCRIPT

cfg = [];
cfg.funparameter = 'pow';
cfg.funcolormap = 'hot';
cfg.funcolorlim = 'zeromax';
% cfg.location = 'max';
if ~isempty(strfind(filename{i}, 'dip05'))
  cfg.location = dip05;
elseif ~isempty(strfind(filename{i}, 'dip06'))
  cfg.location = dip06;
elseif ~isempty(strfind(filename{i}, 'dip07'))
  cfg.location = dip07;
elseif ~isempty(strfind(filename{i}, 'dip08'))
  cfg.location = dip08;
end
ft_sourceplot(cfg, source_lcmv_relative_int);

pngfile = sprintf('%s_lcmv.png', filename{i});
print(pngfile, '-dpng')

matfile = sprintf('%s_lcmv.mat', filename{i});
save(matfile, 'source_lcmv_relative', 'source_lcmv_relative_int');
