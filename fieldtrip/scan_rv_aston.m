% Residual Variance scan

fprintf(fid, '%d, %d, ', amp(i), dip(j));
       
%%

cfg  = [];
cfg.headmodel = headmodel;
cfg.grid = sourcemodel;
cfg.method = 'rv';
cfg.channel = 'megplanar';
% cfg.latency = 0.060;
source_rv_pow = ft_sourceanalysis(cfg, timelock_active);

%%

cfg = [];
cfg.parameter = 'rv'; % converted on the fly to new source structure
cfg.operation = '-log10(x1)';
source_rv = ft_math(cfg, source_rv_pow);

%%

cfg = [];
cfg.funparameter = 'rv';
cfg.location = 'max';
ft_sourceplot(cfg, source_rv);

%%
[~, ind] = max(source_rv.rv);

error_pos = 1e3*norm(source_rv.pos(ind, :)-pos(j, :));

[u, s, v] = svd(source_rv_pow.avg.mom{ind});

amplitude = sqrt(source_rv_pow.avg.pow(ind)) * 10^9;

est_ori = u(:,1)'; % orientation

error_ori = atan2(norm(cross(est_ori, ori(j, :))),dot(est_ori, ori(j, :)));

error_ori = min(abs(error_ori), abs(error_ori-pi));

fprintf(fid, '%f,', [error_ori, error_pos, amplitude, est_ori, source_rv.pos(ind, :), nan]);

if ~isempty(strfind(suffix, 'tsss'))
    fprintf(fid, 'TRUE, ');
else
    fprintf(fid, 'FALSE, ');
end


fprintf(fid, 'ft_rv\n');



%{

cfg.downsample = 2; % for memory reasons
cfg.parameter = 'rv';
source_rv_int = ft_sourceinterpolate(cfg, source_rv, mri);

%%
cfg = [];
% PLEASE NOTE THAT i IS COMING FROM THE CALLING SCRIPT

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

%}