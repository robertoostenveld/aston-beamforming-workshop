% DICS scan

fprintf(fid, '%d, %d, ', amp(i), dip(j));
       
%%

cfg  = [];
cfg.headmodel = headmodel;
cfg.grid = sourcemodel;
cfg.method = 'dics';
cfg.channel = 'megplanar';
cfg.frequency = 20;
cfg.dics.keepcsd = 'yes';
cfg.dics.lambda = '0.01%';
source_dics_active   = ft_sourceanalysis(cfg, freq_active);
source_dics_baseline = ft_sourceanalysis(cfg, freq_baseline);

%%

% avoid rounding error
source_dics_active.freq = round(source_dics_active.freq);
source_dics_baseline.freq = round(source_dics_baseline.freq);

cfg = [];
cfg.parameter = 'pow';
cfg.operation = 'log10(x1/x2)';
source_dics_relative = ft_math(cfg, source_dics_active, source_dics_baseline);

%%
[~, ind] = max(source_dics_relative.pow);

error_pos = 1e3*norm(source_dics_relative.pos(ind, :)-pos(j, :));

[u, s, v] = svd(source_dics_active.avg.csd{ind});

amplitude = sqrt(source_dics_active.avg.pow(ind)) * 10^9;

est_ori = real(u(:,1)'); % orientation

error_ori = atan2(norm(cross(est_ori, ori(j, :))),dot(est_ori, ori(j, :)));

error_ori = min(abs(error_ori), abs(error_ori-pi));

fprintf(fid, '%f,', [error_ori, error_pos, amplitude, est_ori, source_dics_relative.pos(ind, :), nan]);

if ~isempty(strfind(suffix, 'tsss'))
    fprintf(fid, 'TRUE, ');
else
    fprintf(fid, 'FALSE, ');
end


fprintf(fid, 'ft_dics\n');

%{
cfg = [];
cfg.funparameter = 'pow';
cfg.location = 'max';
ft_sourceplot(cfg, source_dics_relative);

%%

cfg = [];
cfg.downsample = 2; % for memory reasons
cfg.parameter = 'pow';
source_dics_relative_int = ft_sourceinterpolate(cfg, source_dics_relative, mri);

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
ft_sourceplot(cfg, source_dics_relative_int);

pngfile = sprintf('%s_dics.png', filename{i});
print(pngfile, '-dpng')

matfile = sprintf('%s_dics.mat', filename{i});
save(matfile, 'source_dics_relative', 'source_dics_relative_int');
%}
