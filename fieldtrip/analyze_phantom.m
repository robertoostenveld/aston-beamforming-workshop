cd C:\home\Data\Aston\phantom

%%
mri = ft_read_mri('CT/Elekta_Vectorview_phantom_ct.nii', 'dataformat', 'nifti_spm');
mri.anatomy = double(mri.anatomy);
mri.anatomy(mri.anatomy==255) = 0;

trans = load('CT/Elekta_Vectorview_phantom_ct.trans');

mri.transform = trans * mri.transform; % FIXME this is not correct
mri = ft_convert_units(mri, 'm');
%%
amp = {'20nAm', '100nAm', '200nAm', '1000nAm'};
dip = {'dip05', 'dip06', 'dip07', 'dip08'};

methods = {'dip_fit', 'scan_rv', 'scan_lcmv', 'scan_dics'};

table = zeros(numel(dip), numel(amp));

% filename = {
%   '1000nAm/dip05_1000nAm_sss.fif'
%   '1000nAm/dip06_1000nAm_sss.fif'
%   '1000nAm/dip07_1000nAm_sss.fif'
%   '1000nAm/dip08_1000nAm_sss.fif'
%   '100nAm/dip05_100nAm_sss.fif'
%   '100nAm/dip06_100nAm_sss.fif'
%   '100nAm/dip07_100nAm_sss.fif'
%   '100nAm/dip08_100nAm_sss.fif'
%   '200nAm/dip05_200nAm_sss.fif'
%   '200nAm/dip06_200nAm_sss.fif'
%   '200nAm/dip07_200nAm_sss.fif'
%   '200nAm/dip08_200nAm_sss.fif'
%   '20nAm/dip05_20nAm_sss.fif'
%   '20nAm/dip06_20nAm_sss.fif'
%   '20nAm/dip07_20nAm_sss.fif'
%   '20nAm/dip08_20nAm_sss.fif'
%   };


% filename = {
%   '1000nAm/dip05_1000nAm.fif'
%   '1000nAm/dip06_1000nAm.fif'
%   '1000nAm/dip07_1000nAm.fif'
%   '1000nAm/dip08_1000nAm.fif'
%   '100nAm/dip05_100nAm.fif'
%   '100nAm/dip06_100nAm.fif'
%   '100nAm/dip07_100nAm.fif'
%   '100nAm/dip08_100nAm.fif'
%   '200nAm/dip05_200nAm.fif'
%   '200nAm/dip06_200nAm.fif'
%   '200nAm/dip07_200nAm.fif'
%   '200nAm/dip08_200nAm.fif'
%   '20nAm/dip05_20nAm.fif'
%   '20nAm/dip06_20nAm.fif'
%   '20nAm/dip07_20nAm.fif'
%   '20nAm/dip08_20nAm.fif'
%   };
%%


%%

dip05 = [37.2 0 52.0]/1000; % 64
dip06 = [27.5 0 46.4]/1000; % 54
dip07 = [15.8 0 41.0]/1000; % 44
dip08 = [7.9  0 33.0]/1000; % 35

%%

headmodel = [];
headmodel.o = [0 0 0];
headmodel.r = 0.10;

%%
pos = [];

for i = 1:numel(filename)
  
  close all
  
  %%
  
  cfg = [];
  cfg.dataset = filename{i};
  cfg.channel = {'megplanar', '-MEG1133', '-MEG1323', '-MEG2233'};
  cfg.trialfun = 'ft_trialfun_general';
  cfg.trialdef.prestim = 0.100;
  cfg.trialdef.poststim = 0.25;
  cfg.trialdef.eventtype = 'STI201';
  cfg.lpfilter = 'yes';
  cfg.lpfreq = 40;
  cfg.demean = 'yes';
  cfg.coilaccuracy = 2;
  cfg = ft_definetrial(cfg);
  data = ft_preprocessing(cfg);
  %%
  cfg = [];
  cfg.method = 'pca';
  cfg.updatesens = 'no';
  cfg.channel = 'megplanar';
  comp = ft_componentanalysis(cfg, data);
  
  cfg = [];
  cfg.updatesens = 'no';
  cfg.component = comp.label(51:end);
  data_fix = ft_rejectcomponent(cfg, comp);
  %%
  cfg = [];
  cfg.toilim = [0 0.1];
  active = ft_redefinetrial(cfg, data_fix);
  
  cfg = [];
  cfg.toilim = [-0.1 0];
  baseline = ft_redefinetrial(cfg, data_fix);
  
  %%
  
  cfg = [];
  cfg.channel = 'meggrad';%'megmag';%'
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
  %ft_multiplotER(cfg, timelock_active);
  
  %%
  
  cfg = [];
  cfg.headmodel = headmodel;
  cfg.grid.resolution = 0.0035;
  sourcemodel = ft_prepare_sourcemodel(cfg, timelock_active);
  
  %%
  cfg=[];
  cfg.headmodel = headmodel;
  cfg.grid.resolution = 20*1e-3;
  cfg.gridsearch = 'yes';
  cfg.channel = 'megplanar';

  cfg.latency = [0.03 0.05];    
  
  dip_source  = ft_dipolefitting(cfg, timelock_active);
     
  %%
  for j = 1:numel(dip)     
          if ~isempty(findstr(dip{j}, filename{i}))
              break;
          end
  end
  
  for k = 1:numel(amp)
      if ~isempty(findstr(amp{k}, filename{i}))
          break;
      end
  end
  %%
  pos = [pos ; dip_source.dip.pos eval(dip{j})];
  
  table(j, k) = 1e3*norm(dip_source.dip.pos-eval(dip{j}));
  %%
  % the remainder is in separate scripts, which generate figures (on disk) and data files
  
  scan_rv
  scan_lcmv
  scan_dics
  
end
%%
save table.mat table dip amp methods