%%
cd /Volumes/ASTON/data/phantom

%%
mri = ft_read_mri('CT/Elekta_Vectorview_phantom_ct.nii');
mri.anatomy = double(mri.anatomy);
mri.anatomy(mri.anatomy==255) = 0;

trans = [
  0.999982 2.86548e-08 -0.0060473 0.00887158
  3.67397e-05 0.999982 0.00606764 0.00284997
  0.00604718 -0.00606776 0.999963 -0.0131477
  0 0 0 1
  ];

mri.transform = trans * mri.transform; % FIXME this is not correct
mri = ft_convert_units(mri, 'm');

filename = {
  '1000nAm/dip05_1000nAm_sss.fif'
  '1000nAm/dip06_1000nAm_sss.fif'
  '1000nAm/dip07_1000nAm_sss.fif'
  '1000nAm/dip08_1000nAm_sss.fif'
  '100nAm/dip05_100nAm_sss.fif'
  '100nAm/dip06_100nAm_sss.fif'
  '100nAm/dip07_100nAm_sss.fif'
  '100nAm/dip08_100nAm_sss.fif'
  '200nAm/dip05_200nAm_sss.fif'
  '200nAm/dip06_200nAm_sss.fif'
  '200nAm/dip07_200nAm_sss.fif'
  '200nAm/dip08_200nAm_sss.fif'
  '20nAm/dip05_20nAm_sss.fif'
  '20nAm/dip06_20nAm_sss.fif'
  '20nAm/dip07_20nAm_sss.fif'
  '20nAm/dip08_20nAm_sss.fif'
  };



%%
dip05 = [32.5 0 56.3]/1000; % 65
dip06 = [27.5 0 47.6]/1000; % 55
dip07 = [22.5 0 39.0]/1000; % 45
dip08 = [17.5 0 30.3]/1000; % 35


%%
headmodel = [];
headmodel.o = [0 0 0];
headmodel.r = 0.10;

%%


i = 1;
cfg = [];
cfg.dataset = filename{i};
cfg.trialfun = 'ft_trialfun_general';
cfg.trialdef.prestim = 0.100;
cfg.trialdef.poststim = 0.25;
cfg.trialdef.eventtype = 'STI201';
cfg.lpfilter = 'yes';
cfg.lpfreq = 40;
cfg.demean = 'yes';
cfg = ft_definetrial(cfg);
data = ft_preprocessing(cfg);

cfg = [];
cfg.toilim = [0 0.100];
active = ft_redefinetrial(cfg, data)

cfg = [];
cfg.toilim = [0.150 0.250];
baseline = ft_redefinetrial(cfg, data)


%%
cfg = [];
cfg.channel = 'megmag';
cfg.viewmode = 'vertical';
ft_databrowser(cfg, data);
% ft_databrowser(cfg, active);
% ft_databrowser(cfg, baseline);

%%

cfg = [];
cfg.covariance = 'yes';
timelock_active = ft_timelockanalysis(cfg, active);
timelock_baseline = ft_timelockanalysis(cfg, baseline);

%%

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.output = 'powandcsd';
cfg.foilim = [0 40];
freq_active = ft_freqanalysis(cfg, active);
freq_baseline = ft_freqanalysis(cfg, baseline);

%%

cfg = [];
cfg.layout = 'neuromag306all.lay';
ft_multiplotER(cfg, timelock_active);

%%

cfg = [];
cfg.headmodel = headmodel;
cfg.grid.resolution = 0.01;
sourcemodel = ft_prepare_sourcemodel(cfg, timelock);

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
cfg.funparameter = 'rv';
cfg.location = 'min';
ft_sourceplot(cfg, source_rv);


%%

cfg = [];
cfg.downsample = 2; % for memory reasons
cfg.parameter = 'rv';
source_rv_int = ft_sourceinterpolate(cfg, source_rv, mri);

%%

cfg = [];
cfg.funparameter = 'rv';
cfg.location = 'min';
cfg.location = dip05;
ft_sourceplot(cfg, source_rv_int);

%%


cfg  = [];
cfg.headmodel = headmodel;
cfg.grid = sourcemodel;
cfg.method = 'lcmv';
cfg.channel = 'megplanar';
cfg.latency = 0.060;
source_rv = ft_sourceanalysis(cfg, active_timelock);

