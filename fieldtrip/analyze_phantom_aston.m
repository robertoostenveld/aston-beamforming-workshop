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

outfile = 'ft_results_nomov.csv';

fid = fopen(outfile, 'w+');
fprintf(fid, ['actual_dipole_amplitude, dipole_index, error_ori, error_pos, estimated_dipole_amplitude, estimated_ori_x, estimated_ori_y, '...
    'estimated_ori_z, estimated_pos_x, estimated_pos_y, estimated_pos_z, gof, maxfilter, method\n']);


amp = [20 200 1000];
dip = 5:12;

methods = {'dip_fit', 'scan_rv', 'scan_lcmv', 'scan_dics'};

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
ori = [
    0.813310537	0	-0.581829846
    0.860261416	0	-0.50985321
    0.93311078	0	-0.359589032
    0.972520896	0	-0.232815608
    0	0.358140538	0.93366769
    0	0.435318817	0.900276362
    0	0.580161604	0.814501389
    0	0.681582014	0.731741729
    
    ];
%%

headmodel = [];
headmodel.o = [0 0 0];
headmodel.r = 0.10;

%%
suffix = '_tsss';%'_movement_tsss_mc';%
%%
outpos = [];

for i = 1:length(amp)
    for j = 1:length(dip)
        
        fprintf(fid, '%d, %d, ', amp(i), dip(j));
        
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
        cfg.foilim = [19 21];
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
        cfg.model  = 'regional';
        
        dip_source  = ft_dipolefitting(cfg, timelock_active);
        
        [u, s, v] = svd(dip_source.dip.mom);
        
        est_ori = u(:,1)'; % orientation
        
        error_ori = atan2(norm(cross(est_ori, ori(j, :))),dot(est_ori, ori(j, :)));
        
        error_ori = min(abs(error_ori), abs(error_ori-pi));
        
        signal = s(1,1)*v(:,1)';
        
        % CHECK Might not be exactly 20Hz
        model = cos(20*2*pi*dip_source.time) + 1i*sin(20*2*pi*dip_source.time);
        
        coeff = signal/model;
        amplitude = 2*abs(coeff); % in Am, zero to peak
        amplitude = 2 * amplitude * 10^9; % peak to peak, in nAm
        %%
        error_pos = 1e3*norm(dip_source.dip.pos-pos(j, :));
        
        gof = 100*(1-var(reshape(dip_source.Vdata-dip_source.Vmodel, 1, []))/ var(reshape(dip_source.Vdata, 1, [])));
        
        fprintf(fid, '%f,', [error_ori, error_pos, amplitude, est_ori, dip_source.dip.pos, gof]);
        
        if ~isempty(strfind(suffix, 'tsss'))
            fprintf(fid, 'TRUE, ');
        else
            fprintf(fid, 'FALSE, ');
        end
        
        
        fprintf(fid, 'ft_dip_fit\n');
        %%
        % the remainder is in separate scripts, which generate figures (on disk) and data files
        
        scan_rv_aston
        scan_lcmv_aston
        scan_dics_aston
        
    end
end
%%
fclose(fid);