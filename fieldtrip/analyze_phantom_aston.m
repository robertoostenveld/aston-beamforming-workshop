cd C:\home\Data\Aston\newphantom

ft_defaults;
%%
% mri = ft_read_mri('CT/Elekta_Vectorview_phantom_ct.nii', 'dataformat', 'nifti_spm');
% mri.anatomy = double(mri.anatomy);
% mri.anatomy(mri.anatomy==255) = 0;
% 
% trans = load('CT/Elekta_Vectorview_phantom_ct.trans');
% 
% mri.transform = trans * mri.transform; % FIXME this is not correct
% mri = ft_convert_units(mri, 'm');
%%
amp = 200;%[20 200 1000];
dip = 5:12;

methods = {'dip_fit', 'scan_rv', 'scan_lcmv', 'scan_dics'};

table = zeros(numel(dip), numel(amp));

%%
pos = [
    37.2 0 52.0;
    27.5 0 46.4;
    15.8 0 41.0;
    7.9  0 33.0;
    0 -59.7 22.9;
    0 -48.6 23.5;
    0 -35.8 25.5;
    0 -24.8 23.1;
    ]/1000;

%%

headmodel = [];
headmodel.o = [0 0 0];
headmodel.r = 0.10;

%%
suffix = '_movement_tsss_mc';%'_tsss';
%%
outpos = [];

for i = 1:length(amp)
    for j = 1:length(dip)
        
        close all
        
        %%
        dataset = sprintf('.\\%dnAm\\Amp%d_Dip%d_IASoff%s.fif',amp(i), amp(i), dip(j), suffix);
        
        hdr = ft_read_header(dataset);
        
        hdr.chantype{strmatch('SYS201', hdr.label)} = 'other trigger';
        
        ev = ft_read_event(dataset, 'header', hdr);
        
        % ev = read_trigger(dataset, 'chanindx', strmatch('SYS201', hdr.label), 'detectflank', 'up');
        
        cfg = [];
        cfg.dataset = dataset;
        cfg.event = ev;
        cfg.channel = {'megplanar', '-MEG1133', '-MEG1323'};        
        cfg.trialfun = 'ft_trialfun_general';
        cfg.trialdef.prestim = 0.100;
        cfg.trialdef.poststim = 0.25;
        cfg.trialdef.eventtype = 'SYS201';
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 40;
        cfg.demean = 'yes';
        %cfg.coilaccuracy = 2;
        cfg = ft_definetrial(cfg);
        data = ft_preprocessing(cfg);
        %%
%         cfg = [];
%         cfg.method = 'pca';
%         cfg.updatesens = 'no';
%         cfg.channel = 'megplanar';
%         comp = ft_componentanalysis(cfg, data);
%         
%         cfg = [];
%         cfg.updatesens = 'no';
%         cfg.component = comp.label(51:end);
%         data_fix = ft_rejectcomponent(cfg, comp);
        %%
        cfg = [];
        cfg.toilim = [0 0.1];
        active = ft_redefinetrial(cfg, data);
        
        cfg = [];
        cfg.toilim = [-0.1 0];
        baseline = ft_redefinetrial(cfg, data);
        
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
        %cfg.channel = 'megmag';
        cfg.layout = 'neuromag306all.lay';
        ft_multiplotER(cfg, timelock_active);
        
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
        outpos = [outpos; dip_source.dip.pos pos(j, :)];
        
        table(j, i) = 1e3*norm(dip_source.dip.pos-pos(j, :));
        %%
        % the remainder is in separate scripts, which generate figures (on disk) and data files
        
        scan_rv_aston
        scan_lcmv_aston
        scan_dics_aston
        
    end
end
%%
save table.mat table dip amp methods