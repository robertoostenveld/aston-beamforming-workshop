% Residual Variance scan

%%

cfg  = [];
cfg.headmodel = headmodel;
cfg.grid = sourcemodel;
cfg.method = 'rv';
cfg.channel = 'megplanar';
% cfg.latency = 0.060;
source_rv = ft_sourceanalysis(cfg, timelock_active);

%%

cfg = [];
cfg.parameter = 'rv'; % converted on the fly to new source structure
cfg.operation = '-log10(x1)';
source_rv = ft_math(cfg, source_rv);

%%

cfg = [];
cfg.funparameter = 'rv';
cfg.location = 'max';
ft_sourceplot(cfg, source_rv);


%%

cfg = [];
cfg.downsample = 2; % for memory reasons
cfg.parameter = 'rv';
source_rv_int = ft_sourceinterpolate(cfg, source_rv, mri);

%%
% PLEASE NOTE THAT i IS COMING FROM THE CALLING script

cfg = [];
cfg.funparameter = 'rv';
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
ft_sourceplot(cfg, source_rv_int);

pngfile = sprintf('%s_rv.png', filename{i});
print(pngfile, '-dpng')

matfile = sprintf('%s_rv.mat', filename{i});
save(matfile, 'source_rv', 'source_rv_int');

