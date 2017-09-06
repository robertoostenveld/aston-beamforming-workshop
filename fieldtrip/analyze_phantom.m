cd C:\home\Data\Aston\phantom

%%
mri = ft_read_mri('CT/Elekta_Vectorview_phantom_ct.nii', 'dataformat', 'nifti_spm');
mri.anatomy = double(mri.anatomy);
mri.anatomy(mri.anatomy==255) = 0;

trans = load('CT/Elekta_Vectorview_phantom_ct.trans');

mri.transform = trans * mri.transform; % FIXME this is not correct
mri = ft_convert_units(mri, 'm');
%%

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

dip05 = [37.2 0 52  ]/1000; % 64
dip06 = [27.5 0 46.4]/1000; % 54
dip07 = [15.8 0 41.0]/1000; % 44
dip08 = [7.9  0 30.3]/1000; % 35

%%

headmodel = [];
headmodel.o = [0 0 0];
headmodel.r = 0.10;

%%

for i = 1:numel(filename)
  
  close all
  
  %%
  
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
  %%
  cfg = [];
  cfg.toilim = [0 0.1];
  active = ft_redefinetrial(cfg, data);
  
  cfg = [];
  cfg.toilim = [-0.1 0];
  baseline = ft_redefinetrial(cfg, data);
  
  %%
  
  cfg = [];
  cfg.channel = 'meggrad';%'megmag';
  cfg.viewmode = 'vertical';
  % ft_databrowser(cfg, data);
  % ft_databrowser(cfg, active);
  % ft_databrowser(cfg, baseline);
  
  %%
  
  cfg = [];
  cfg.covariance = 'yes';
  cfg.covariancewindow = 'all';
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
  cfg.channel = 'megmag';
  cfg.layout = 'neuromag306all.lay';
  ft_multiplotER(cfg, timelock_active);
  
  %%
  
  cfg = [];
  cfg.headmodel = headmodel;
  cfg.grid.resolution = 0.01;
  sourcemodel = ft_prepare_sourcemodel(cfg, timelock_active);
  
  %%
  % the remainder is in separate scripts, which generate figures (on disk) and data files
  
  scan_rv
  scan_lcmv
  scan_dics
  
end
