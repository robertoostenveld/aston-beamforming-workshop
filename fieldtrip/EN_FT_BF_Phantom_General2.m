%% Fieldtrip beamforming pipeline for phantom data
% Date: 21/02/2017
% Author: Amit @ Elekta Neuromag
% Details: Code for all phantom data
%% Add fieldtrip in path
clear all; clc; refresh; close all;
restoredefaultpath
add_ft='yes'; add_spm='no';
addpath('~//Dropbox//MATLAB//')
homedir=EN_add_toolboxes(add_ft, add_spm);
cd('~//spm_temp//')
tic;

%% Set parameters for the analysis
allsite={'Aalto', 'Aston', 'Aston_mov', 'Bari' 'Biomag'}; % Aalto Aston
datacat='phantom';
chancat=[1 2 3]; 
SSS=[0 1];
apply_lcmv = 'yes';
apply_dics = 'no';
dics_freq = [5.0, 35.0];
reg_form = 'SNR^4/500';

%%  Run the loop                          
for sitenum=1:5
    site=allsite{1,sitenum};
    if strcmp(site, 'Aalto')
        maxf=[1 2]; amps=[2,4,5,7]; dipoles= 5:8;
    elseif strcmp(site, 'Aston')
        maxf=[1 2]; amps=[2 5 7]; dipoles= 5:12;
    elseif strcmp(site, 'Bari')
        maxf=[1 2]; amps=3; dipoles = 5:12;
    elseif strcmp(site, 'Biomag')
        maxf=[1 2]; amps=6; dipoles = 1:12;
    elseif strcmp(site, 'Aston_mov')
        maxf=[1 4]; amps=5; dipoles=5:12; phantom_moving='yes';
    end
    
    [actual_diploc, channelslist]=EN_load_plantom_actual_dip_chanlist(site);
%% Main loop
    for meg = chancat(1)
        for pp = maxf
            for S = amps
                for dip = dipoles
                    
% Set parameters for the analysis
    par.site            = site;
    par.Movephantom     = 'no';
    par.prep            = {'', '_sss', '_tsss', '_tsss_mc', '_cxsss', '_nosss'};
    par.sub             = {'10', '20', '25', '100', '200', '500', '1000'};
    par.dipoles         = 1:32;
    par.meg             = {'all', 'mag', 'grad'};
    par.stimraise       = 0;
    par.apply_dipfit    = 'yes';
    par.apply_lcmv      = 'yes';
    par.apply_dics      = 'no';
    par.reg             = 5; % default
    par.dics_freq       = [5.0, 35.0];
    par.reg_compare     = '';
    par.visual          = 'no';
    par.powspect        = '';
    par.browse          = ''; 
    par.more_plots      = '';
    par.numtry          = 1;
    par.gridres         = 0.005;
    par.bpfreq          = [2 95];
    par.resultplot      = '';
    par.resultplotsave  = '';

%% Set data directory and other parameters
    [workdir, megfile, out_path]=EN_set_directory2(site, datacat, par.sub, S, dip, par.prep, pp);
    out_dir=[out_path  'FieldTrip//'];
    phantom_ct_file='/neuro/mri/phantom-ct-1/sets/fantomi_1362_01-jne-100423-5.nii';
    [mdir, mfname,~]=fileparts(megfile);
    mkdir(out_dir)
    
    if isequal(par.meg{1,meg}, 'all')
        ignorech = [];
    elseif isequal(par.meg{1,meg}, 'mag')
        ignorech = channelslist.grad;
    elseif isequal(par.meg{1,meg}, 'grad')
        ignorech = channelslist.mag;
    end
    % Aalto Phatom data
    if isequal(par.site,'Aalto') 
        if isequal(par.prep{1,pp}, '')
            bads       = {'MEG0111', 'MEG0511', 'MEG2231', 'MEG2233', 'MEG2422', 'MEG1842', 'MEG0433'}; 
        else
            bads       = {};
        end
        badch           = strcat('-', [ignorech, bads]);
        megchan         = ['MEG*', badch];
        stimchan        = 'STI201';
        par.trig_min_gap= 0.1;
        par.trialwin    = [-0.1 0.1];
        par.ctrlwin     = [-0.1 0.0];
        par.actiwin     = [0.0 0.1];

    % Aston Phatom data
    elseif isequal(par.site,'Aston') || isequal(par.site,'Aston_mov')
        if isequal(par.prep{1,pp}, '')
            bads       = {'MEG0613', 'MEG1032', 'MEG1133', 'MEG1323'}; 
        else
            bads       = {};
        end
        badch            = strcat('-', [ignorech, bads]);
        megchan          = ['MEG*', badch];
        stimchan         = 'SYS201';
        par.trig_min_gap = 0.5;
        par.trialwin     = [-0.5 0.5];
        par.ctrlwin      = [-0.5 0.0];
        par.actiwin      = [0.0 0.5];
        if isequal(par.site, 'Aston_mov')
            par.stimraise   = 3840;
            par.Movephantom = 'yes';
        end

    % Biomag Phatom data
    elseif isequal(par.site,'Biomag') 
        if isequal(par.prep{1,pp}, '')
            bads       = {'MEG0342',  'MEG0542', 'MEG1013'}; 
        else  
            bads       = {};
        end
        badch            = strcat('-', [ignorech, bads]);
        stimchan         = 'SYS201';
        megchan          = ['MEG*', badch];
        par.trig_min_gap = 0.1;
        par.trialwin     = [-0.1 0.1];
        par.ctrlwin      = [-0.1 0.0];
        par.actiwin      = [0.0 0.1];
        
    % Bari Phatom data
    elseif isequal(par.site,'Bari')
        if isequal(par.prep{1,pp}, '') || isequal(par.prep{1,pp}, '_nosss')
            bads       = {'MEG0943', 'MEG0631', 'MEG0222', 'MEG1522', 'MEG1512', 'MEG1432', 'MEG1113'}; 
        else  
            bads       = {};
        end
        badch            = strcat('-', [ignorech, bads]);
        megchan          = ['MEG*', badch];
        stimchan         = 'SYS201';
        par.trig_min_gap = 0.1;
        par.trialwin     = [-0.1 0.1];
        par.ctrlwin      = [-0.1 0.0];
        par.actiwin      = [0.0 0.1];
    end

%% Find trigger categories && label them
    keyset={['Dipole-' strcat(num2str(dip))]};
    valueset=[dip];
    evdict=containers.Map(keyset, valueset);
%% Browse raw data
if isequal(par.browse,'yes')
    cfg          = [];
    cfg.dataset  = megfile;
    cfg.channel  = megchan;
    cfg.viewmode = 'butterfly';
    cfg.blocksize= 15;
    cfg.checkmaxfilter      = 0;
    cfg.demean   = 'yes';
    ft_databrowser(cfg);
end       
%% Prepare data 
% Define trials
    cfg                     = [];                   
    cfg.dataset             = megfile;
    cfg.channel             = megchan;
    cfg.checkmaxfilter      = 0;
    cfg.trialdef.eventtype  = stimchan;
    cfg.trialdef.eventvalue = dip + par.stimraise;                     
    cfg.trialdef.prestim    = abs(par.trialwin(1)) + 0.01;                     
    cfg.trialdef.poststim   = par.trialwin(2) + 0.01; 
    cfg.trialfun            = 'EN_phantom_trialfun';
    cfg.trig_min_gap        = par.trig_min_gap;
    cfg.minlength = cfg.trialdef.prestim + cfg.trialdef.poststim;
    cfg = ft_definetrial(cfg);
    if isequal(par.Movephantom, 'yes')
        cfg.trl(:,4)=cfg.trl(:,4)-par.stimraise;
    end
% Leave the first and last trigger
        cfg.trl(1:1,:)  =[];
        cfg.trl(end:end,:)=[];

% Data preprocessing  
    cfg.demean          = 'yes';
    cfg.baselinewindow  = [par.trialwin(1) 0.0];
    cfg.detrend         = 'yes';
    cfg.bpfilter        = 'yes'; 
    cfg.bpfilttype      = 'but';
    cfg.bpfreq          = par.bpfreq;
    cfg.coilaccuracy    = 1;
    cfg.checkmaxfilter  = 0;
    data = ft_preprocessing(cfg);

    cfg.bsfilter='yes';
    cfg.bsfreq = [49.5 50.5];
    data = ft_preprocessing(cfg, data);
    
%% Resegment for the actual trial window
    cfg2=[];
    cfg2.toilim = par.trialwin;
    data = ft_redefinetrial(cfg2, data);

%% Visualize events triggers && number of trials per category
    EN_ft_plot_events(figure, cfg, keyset, valueset)

%% Interactive data browser 
if isequal(par.more_plots, 'yes')
    cfg.blocksize          = cfg.trialdef.prestim + cfg.trialdef.poststim;
    cfg.continuous         = 'no';
    cfg.viewmode           = 'butterfly'; %''
    ft_databrowser(cfg, data);
end

%% Visual bad trial detection and rejection
if isequal(par.visual, 'yes')  
    cfg          = [];
    cfg.method   = 'summary';
    cfg.metric   = 'min';
    data_summary = ft_rejectvisual(cfg,data);
else
    data_summary = data;
end
    clear('data')
%% Trial averaging and computation of covariance matrix
    cfg = [];
    cfg.covariance='yes';
    cfg.covariancewindow = par.trialwin; %'all';
    cfg.vartrllength = 2;
    timelock = ft_timelockanalysis(cfg,data_summary);
    
    if isequal(par.more_plots, 'yes')
        figure
        subplot_tight(3,4,1:4,0.05), plot(timelock.time, timelock.avg);
        xlim([timelock.time(1) timelock.time(end)]), title('Averaged timelocked response')   
        subplot_tight(3,4,[5 6 9 10],0.05), cfg.layout = timelock.grad; %'neuromag306all.lay';
        ft_multiplotER(cfg, timelock); title('Evoked response topoplot')
        subplot_tight(3,4,[7 8 11 12],0.05),imagesc(timelock.cov)
        title('Covariance plot for timelock')
    end
%% Redefine pre and post data and take average
    cfg = [];
    cfg.toilim = par.actiwin;
    datapst = ft_redefinetrial(cfg, data_summary);
    cfg = [];
    cfg.covariance='yes';
    avgpst = ft_timelockanalysis(cfg,datapst);

    cfg = [];            
    cfg.toilim = par.ctrlwin;
    datapre = ft_redefinetrial(cfg, data_summary);
    cfg = [];
    cfg.covariance='yes';
    avgpre = ft_timelockanalysis(cfg,datapre);
    
    ntrials = length(datapst.trial);
    SNR=snr(avgpst.avg',avgpre.avg');

    par.reg = eval(reg_form);
    
%% Plot redefined data and their covariance
    if isequal(par.more_plots, 'yes')
        figure %('units', 'normalized', 'outerposition', [0 0 1 1])
        subplot_tight(2,2,1,0.05); plot(avgpre.time, avgpre.avg); xlim([avgpre.time(1) avgpre.time(end)]) 
        subplot_tight(2,2,2,0.05); plot(avgpst.time, avgpst.avg); xlim([avgpst.time(1) avgpst.time(end)])
        subplot_tight(2,2,3,0.05); imagesc(avgpre.cov)
        subplot_tight(2,2,4,0.05); imagesc(avgpst.cov)
    end
    clear datapre datapst data_summary hdr
%% Define anatomical volume or headmodel
    headmodel=[];
    headmodel.o       = [0,0,0];  % in cm
    headmodel.r       = .070;      % should not matter
    headmodel= ft_convert_units(headmodel, 'm'); % prepare head volume
    % ft_plot_mesh(ft_prepare_mesh([], headmodel)); camlight; rotate3d

%% Dipole fitting
if isequal(par.apply_dipfit, 'yes')
    clear max; [maxx, indd]=max(abs(avgpst.avg(:))); % taking max of abs values 
    [ch_peak, t_peak]=ind2sub(size(avgpst.avg), indd);
    timepoint=avgpst.time(t_peak);  % peak timepoint for dipfit

    cfg=[];
    cfg.headmodel   = headmodel;
    cfg.grid.xgrid  = -0.065:par.gridres:0.065;
    cfg.grid.ygrid  = -0.065:par.gridres:0.065;
    cfg.grid.zgrid  = -0.000:par.gridres:0.065;
    cfg.inwardshift = 0.0025;
    cfg.gridsearch  = 'yes';
    cfg.channel     = megchan;
    cfg.latency     = [timepoint,  timepoint]; 
    cfg.numdipoles  = 1;
    dip_source  = ft_dipolefitting(cfg, avgpst);

    hspot = mean(dip_source.dip.pos, 1)*1000;
    hval= max(max(dip_source.dip.mom))*1e9;
    difff=sqrt(sum((actual_diploc(dip,:)-hspot).^2));
    
    disp([hspot, hval, difff])

    resultfile=[out_dir 'FT' '-' site '_DipoleFitting_results.csv'];
    fid = fopen(strcat(resultfile), 'a+');
    fprintf(fid, '%s,', mfname);
    fprintf(fid, '%s,', char(par.sub(S)));            
    fprintf(fid, '%f,', [dip, hspot(1),hspot(2),hspot(3), hval, difff]);
    fprintf(fid, '%s,', par.prep{1,pp});
    fprintf(fid, '%s,', 'LCMV');
    fprintf(fid, '\n');
    fclose(fid);
            
    if isequal(par.more_plots, 'yes')
        if ~exist('phantom_ct', 'var')
            phantom_ct=ft_convert_units(ft_read_mri(phantom_ct_file), 'm');
        end
        figure
        hold on
        ft_plot_dipole(dip_source.dip.pos(1,:), mean(dip_source.dip.mom(1:3,:),2), 'color', 'r', 'unit', 'm')
        pos = mean(dip_source.dip.pos, 1);
        ft_plot_slice(phantom_ct.anatomy, 'transform', phantom_ct.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.001)
        ft_plot_slice(phantom_ct.anatomy, 'transform', phantom_ct.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.001)
        ft_plot_slice(phantom_ct.anatomy, 'transform', phantom_ct.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.001)
        ft_plot_crosshair(pos, 'color', [1 1 1]/2);
        title(['Dipole fitting plot for ' mfname '-dip' strcat(num2str(dip))], 'FontSize', 15)
        axis tight
        view(-15,15)
        rotate3d
    end 
end
clear dip_source

%% Apply LCMV beamformer
if isequal(par.apply_lcmv, 'yes')
% Prepare/Compute leadfield
    cfg             = [];
    cfg.grad        = ft_convert_units(timelock.grad, 'm');  
    cfg.headmodel   = headmodel;  
    cfg.grid.xgrid  = -0.065:par.gridres:0.065; % Define rectangular grid
    cfg.grid.ygrid  = -0.065:par.gridres:0.065;
    cfg.grid.zgrid  = -0.000:par.gridres:0.065;  
    cfg.channel     = timelock.label;
    cfg.normalize   = 'yes'; %remove depth bias 
    leadfield       = ft_prepare_leadfield(cfg);
            
    if isequal(par.more_plots, 'yes')
        figure
        ft_plot_headshape(ft_convert_units(ft_read_headshape(megfile), 'm'))
        ft_plot_vol(headmodel, 'facecolor', 'skin'); alpha 0.3; camlight
        ft_plot_mesh(leadfield.pos(leadfield.inside,:), 'vertexcolor', 'skin')
        ft_plot_sens(ft_convert_units(ft_read_sens(megfile), 'm'), 'coilshape', 'circle', 'coilsize', 0.015, 'facecolor', [1 1 1])
        rotate3d 
    end

%% Source analysis
% Creat spatial filter using LCMV beamformer
    cfg                  = [];
    cfg.method           = 'lcmv';
    cfg.grid             = leadfield; % leadfield, which has the grid information
    cfg.headmodel        = headmodel; % volume conduction model (headmodel)
    cfg.lcmv.keepfilter  = 'yes';
    cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
    cfg.lcmv.lambda      = [num2str(par.reg) '%'];
    source_avg           = ft_sourceanalysis(cfg, timelock); % ft_sourceanalysis to create spatial filters
% Apply computed spatial filters on pre and post data     
    cfg                  =[];
    cfg.method           ='lcmv';
    cfg.grid             =leadfield;
    cfg.grid.filter      =source_avg.avg.filter;
    cfg.headmodel        =headmodel;
    sourcepreM1=ft_sourceanalysis(cfg, avgpre); % ft_sourceanalysis to apply spatial filter on pre data
    sourcepstM1=ft_sourceanalysis(cfg, avgpst); % ft_sourceanalysis to apply spatial filter on post data

%% Plot source time series (stc in MNE-p way)
if isequal(par.more_plots, 'yes')
    n_max = 5;
    virt_pow=sourcepreM1.avg.pow;
    virt_pow(isnan(virt_pow)) = 0;
    [descend_val,descend_pos] = sort(virt_pow, 'descend'); % arange in descending order
    nmax_values = descend_val(1:n_max);
    nmax_pos = descend_pos(1:n_max);
    stcpre=cell2mat(sourcepreM1.avg.mom(nmax_pos))';
    virt_pow=sourcepstM1.avg.pow;
    virt_pow(isnan(virt_pow)) = 0;
    [descend_val,descend_pos] = sort(virt_pow, 'descend'); % arange in descending order
    nmax_values = descend_val(1:n_max);
    nmax_pos = descend_pos(1:n_max);
    stcpst=cell2mat(sourcepstM1.avg.mom(nmax_pos))';
    figure, plot([sourcepreM1.time, sourcepstM1.time], [abs(stcpre); abs(stcpst)], 'linewidth', 2)
    title(['STC plot for ' strcat(num2str(n_max)) ' strongest sources'], 'FontSize', 12)
end 
%% Calculate and save NAI (Overlay)
    M1          = sourcepstM1;
    M1.avg.pow  = (sourcepstM1.avg.pow-sourcepreM1.avg.pow)./sourcepreM1.avg.pow;
%% Source interpolate
    cfg              = [];
    cfg.voxelcoord   = 'no';
    cfg.parameter    = 'avg.pow';
    cfg.interpmethod = 'nearest';
    source_int  = ft_sourceinterpolate(cfg, M1, leadfield);
    source_int_mm =  ft_convert_units(source_int, 'mm');
            
%% Write hotspot (dipole) location in a text file
    [~,hind] = max(abs(source_int_mm.pow));
    hval=source_int_mm.pow(hind);
    hspot=source_int_mm.pos(hind,:);
    difff=sqrt(sum((actual_diploc(dip,:)-hspot).^2));
    disp([hspot, hval, difff])

    resultfile=[out_dir 'FT_BF' '-' site '_estimated-source_loc.csv'];
    fid = fopen(strcat(resultfile), 'a+');
    fprintf(fid, '%s,', mfname);
    fprintf(fid, '%s,', char(par.sub(S)));
    fprintf(fid, '%s,', [num2str(dip) '/' char(par.sub(S)) '/' par.prep{1,pp}(2:end)]);
    fprintf(fid, '%f,', [hspot(1),hspot(2),hspot(3), hval, difff, ntrials, size(avgpst.avg,1), SNR, par.reg]);
    fprintf(fid, '%s,', 'LCMV');
    fprintf(fid, '\n');
    fclose(fid);
            
%% Source plot
    if isequal(par.resultplot, 'yes')
        phantom_ct=ft_convert_units(ft_read_mri(phantom_ct_file), 'mm');
        source_int_mm.mask = source_int_mm.pow > max(source_int_mm.pow(:))*0.20; % Set threshold for plotting
        cfg                 = [];
        cfg.method          = 'ortho';
        cfg.funparameter    = 'avg.pow';
        cfg.downsample      = 5;
        cfg.maskparameter   = 'mask';
        cfg.funcolormap     = 'jet';
        cfg.colorbar        = 'yes';
        ft_sourceplot(cfg, M1, phantom_ct);
        camroll(180)
        title(['Source plot for ' mfname '-' par.meg{1,meg}], 'FontSize', 15) %
        if isequal(par.resultplotsave, 'yes')
            set(gcf, 'Position', get(0, 'Screensize'));
            mkdir([out_dir 'FieldTrip_results/'])
            saveas(gcf, [out_dir 'FieldTrip_results//' mfname '-' par.meg{1,meg}], 'tiff')
            close(gcf)
        end
    end


end
    close all;
    clearex par meg pp S dip actual_diploc hspot16 site maxf amps dipoles sitenum ...
            allsite datacat chancat SSS apply_dics apply_lcmv reg_compare dics_freq ...
            phantom_moving reg_form
    toc;
                end
            end
        end
    end
end
datetime
toc
  
