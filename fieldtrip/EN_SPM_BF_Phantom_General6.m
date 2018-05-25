%% SPM beamforming pipeline (Combining LCMV and DICS)
%% Date: 23/11/2017
%% Author: Vladimir @ UCL && Amit @ Elekta Neuromag

%% Add fieldtrip and spm in path
clear all; clc; refresh; close all;
restoredefaultpath
add_ft='yes'; add_spm='yes';
addpath('~//Dropbox//MATLAB//')
homedir=EN_add_toolboxes(add_ft, add_spm);
cd('~//spm_temp//')
tic;

%% Set parameters for the analysis
allsite = {'Aalto', 'Aston', 'Aston_mov', 'Bari', 'Biomag'}; 
datacat = 'phantom';
chancat = [1]; % 1='ALL', 2='MAG', 3='GRAD'
SSS = [0 1];
apply_lcmv = 'yes';
apply_dics = 'no';
reg_compare = 'no';
dics_freq = [5.0, 35.0];
reg_form = 'SNR^4/500';

%%  Run the loop                          
for iii=1:5
    site=allsite{1,iii};

    if strcmp(site, 'Aalto')
        maxf=[1 2]; amps=[1 3 4 6]; dipoles= 5:8;
    elseif strcmp(site, 'Aston')
        maxf=[1 2]; amps=[1 4 6];   dipoles= 5:12;
    elseif strcmp(site, 'Bari')
        maxf=[2];   amps=2;         dipoles =5:12;
    elseif strcmp(site, 'Biomag')
        maxf=[1 2]; amps=5;         dipoles = 5:12;
    elseif strcmp(site, 'Aston_mov')
        maxf=[1 4]; amps=4;         dipoles=5:12;   phantom_moving = 'yes';
    end
    
    if isequal(site, 'Aalto') % For Vectorview
        actual_diploc=load('~/Documents/MATLAB/aalto_phantom.mat'); 
        actual_diploc=actual_diploc.aalto_phantom;
    elseif any(strcmp({'Aston', 'Aston_mov', 'Bari', 'Biomag'}, site)) % For Triux
        actual_diploc=load('~/Documents/MATLAB/biomag_phantom.mat');
        actual_diploc=actual_diploc.biomag_phantom;
    else
    end

%% Main loop
    for meg = chancat
        for pp = maxf
            for sss = SSS(1)
                for sub = amps
                    for dip = dipoles
% Set data directory and other parameters
        par            =[];
        par.subjects   = {'20', '25','100', '200', '500', '1000'};
        par.prep       = {'', '_sss', '_tsss', '_tsss_mc','_cxsss', '_nosss'};
        par.meg        = {'all', 'mag', 'grad'};
        par.visual     = 'no';
        par.powspect   = 'no';
        par.browse     = 'no'; 
        par.more_plots = 'no';
        par.runtry     = 1;
        par.gridres    = 7.0; % in mm
        par.bpfreq     = [2 95];
        par.bsfreq     = [49.1 50.9];
        par.stimraise  = 0;
        
        [workdir, megfile, out_path]=EN_set_directory2(site, datacat, par.subjects, sub, dip, par.prep, pp);
        out_dir=[out_path  'SPM/'];
        mrifile ='/neuro/data/beamformer_data/Beamformer_share/phantom/CT/fantomi_1362_01-jne-100423-5.nii,1';
        load([homedir 'Dropbox/MATLAB/channelslist.mat']);
        if isequal(par.meg{1,meg}, 'all')
            ignorech = [];
        elseif isequal(par.meg{1,meg}, 'mag')
            ignorech = channelslist.grad;
        elseif isequal(par.meg{1,meg}, 'grad')
            ignorech = channelslist.mag;
        end

        if isequal(site, 'Aalto') % For Aalto phantom data  
            if isequal(par.prep{1,pp}, '')
                bads       = {'MEG0111', 'MEG0511', 'MEG2231', 'MEG2233', 'MEG2422', 'MEG1842', 'MEG0433'}; 
            else
                bads       = {};
            end
            badch           = [ignorech, bads];
            par.stimchan    = 'STI201';
            par.stimtype    = [par.stimchan '_up'];
            par.trig_min_gap= 0.1;
            par.trial_win   = [-0.1 0.1];
            par.ctrlwin_ph  = [-0.1 0.0];
            par.actiwin_ph  = [0.0 0.1];
            par.win_oi      = [-0.100 0.000; 0.000 0.100];
            par.badtrs=[];

        elseif isequal(site, 'Aston') % For Aston phantom data     
            if isequal(par.prep{1,pp}, '')
                bads       = {'MEG0613', 'MEG1032', 'MEG1133', 'MEG1323'}; 
            else
                bads       = {};
            end
            badch = [ignorech, bads];
            par.stimchan = 'SYS201';
            par.stimtype   = [par.stimchan '_up'];
            par.trig_min_gap = 0.5;
            par.trial_win = [-0.5 0.5];
            par.actiwin_ph = [0.0 0.5];
            par.ctrlwin_ph  =[-0.5 0.0];
            par.win_oi     = [-0.500 0.000; 0.000 0.500];
            par.badtrs=[];
            
        elseif isequal(site, 'Aston_mov')
            if isequal(par.prep{1,pp}, '')
                bads       = {'MEG0613', 'MEG1032', 'MEG1133', 'MEG1323'}; 
            else
                bads       = {};
            end
            badch = [ignorech, bads];
            par.stimchan = 'SYS201';
            par.stimtype   = [par.stimchan '_up'];
            par.trig_min_gap = 0.5;
            par.trial_win = [-0.5 0.5];
            par.actiwin_ph = [0.0 0.5];
            par.ctrlwin_ph  =[-0.5 0.0];
            par.win_oi     = [-0.500 0.000; 0.000 0.500];
            par.stimraise = 3840;
            par.badtrs=[];
            
        elseif isequal(site,'Biomag') % For Biomag phantom data  
            if isequal(par.prep{1,pp}, '')
                bads       = {'MEG0342',  'MEG0542', 'MEG1013'}; 
            else  
                bads       = {};
            end
            badch          = [ignorech, bads];
            par.stimchan   = 'SYS201';
            par.stimtype   = [par.stimchan '_up'];
            par.trig_min_gap = 0.1;
            par.trial_win  = [-0.1 0.1];
            par.actiwin_ph = [0.0 0.1];
            par.ctrlwin_ph = [-0.1 0.0];
            par.win_oi     = [-0.100 0.000; 0.000 0.100];
            par.badtrs     = [];

        elseif isequal(site,'Bari') % For Bari phantom data                          
            if isequal(par.prep{1,pp}, '')
                bads       = {'MEG0943',  'MEG0222', 'MEG1522', 'MEG1512', 'MEG1432', 'MEG1113', 'MEG0631'}; 
            else  
                bads       = {};
            end
            badch          = [ignorech, bads];
            par.stimchan   = 'SYS201';
            par.stimtype   = [par.stimchan '_up'];
            par.trig_min_gap= 0.1;
            par.trial_win  = [-0.1 0.1];
            par.ctrlwin_ph = [-0.1 0.0];
            par.actiwin_ph = [0.0 0.1];
            par.win_oi     = [-0.100 0.000; 0.000 0.100];
            par.badtrs     = [];%1:15;
        end
        
        [mdir, mfname, ~]=fileparts(megfile);
          
        megchan=channelslist.meg;
        megchan(ismember(megchan, badch))=[]; % Make a channel list without bads + ignorech
        
        trialwin= par.trial_win*1000;  % convert into milisecond
        woi     = par.win_oi*1000;     % convert into milisecond

%% Browse raw in FieldTrip
    if isequal(par.browse, 'yes')
        cfg             = [];
        cfg.dataset     = megfile;
        cfg.channel     = megchan;
        cfg.viewmode    = 'vertical';
        cfg.blocksize   = 1;
        cfg.demean      = 'yes';
        ft_databrowser(cfg);    
    end  
%% Convert fif raw to SPM epoched format 
        S = [];
        S.dataset = megfile;
        S.mode = 'epoched';
        S.channels = megchan;
        S.saveorigheader = 1;
        S.inputformat  ='neuromag_fif'; 
        S.conditionlabels=['dip' num2str(dip)];
        S.trl=EN_make_trialdef_for_SPMbatch5(megfile, par.stimchan,...
                                    dip+par.stimraise, par.trial_win, par.trig_min_gap, 1, 0);
        D = spm_eeg_convert(S);
                
 %% Remove bad trials, if any 
        if par.badtrs
            D = badtrials(D, par.badtrs, 1);
            save(D);

            S=[];
            S.D=D;
            S.prefix='r';
            D=spm_eeg_remove_bad_trials(S);
        end
        
%% Apply bandpass and notch filter
        S       = [];
        S.D     = D;
        S.type  = 'butterworth';
        S.band  = 'bandpass';
        S.freq  = par.bpfreq;
        S.dir   = 'twopass';
        S.prefix= 'f';
        D = spm_eeg_filter(S);

        S       = [];
        S.D     = D;
        S.type  = 'butterworth';
        S.band  = 'stop';
        S.freq  = par.bsfreq;
        S.dir   = 'twopass';
        S.prefix= 'f';
        D = spm_eeg_filter(S);       
        
%% Browse converted raw data
    if isequal(par.browse, 'yes')
        cfg             = [];
        cfg.channel     = megchan;
        cfg.viewmode    = 'vertical';
        cfg.blocksize   = abs(par.trial_win(1)) + par.trial_win(2);
        cfg.demean      = 'yes';
        ft_databrowser(cfg, D.ftraw);
        
        cfg.viewmode    = 'butterfly';
        ft_databrowser(cfg, D.ftraw);
    end  

%%  Power Spectrum analysis and visualization (Additional)
    if isequal(par.powspect, 'yes')
        cfg              = [];
        cfg.output       = 'pow'; 
        cfg.channel      =  megchan;
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.toi          = par.trial_win(1):0.01:par.trial_win(2);
        cfg.foi          = par.bpfreq(1):2:par.bpfreq(2);
        cfg.t_ftimwin    = ones(size(cfg.foi)) * 0.01;
        cfg.trials       = find(D.ftraw.trialinfo(:,1) == dip);
        TFR     = ft_freqanalysis(cfg, D.ftraw);

        cfg = [];
        cfg.baseline     =[par.trial_win(1) 0.0];
        cfg.baselinetype = 'absolute'; 
        cfg.layout       = D.sensors('MEG');
        cfg.showlabels   = 'no';  
        figure 
        cfg.layout       = 'neuromag306planar.lay';
        subplot(1,2,1), ft_multiplotTFR(cfg, TFR);
        cfg.layout       = 'neuromag306mag.lay';
        subplot(1,2,2), ft_multiplotTFR(cfg, TFR);
    end   
    clear TFR

%% Apply TSSS on raw data
    if sss
        disp(channelslist)
        S = [];
        S.D = D;
        S.tsss       = 1;     
        S.t_window   = 1;     
        S.corr_limit = 0.98;  
        S.magscale   = 59.5;     
        S.xspace     = 0;     
        S.Lin        = 8;     
        S.Lout       = 3;     
        S.cond_threshold = 50; 
        S.prefix     = 'sss_'; 
        D = tsss_spm_enm(S);

        D = badchannels(D, ':', 0);save(D);

        if isequal(par.more_plots, 'yes')
            cfg             = [];
            cfg.viewmode    = 'butterfly';
            cfg.blocksize   = 15;
            ft_databrowser(cfg, D.ftraw);     
        else
        end 

    else
        D = montage(D, 'switch', 0);save(D);
    end
    
% Run the batch    
        clear matlabbatch
        if sss
            matlabbatch{1}.spm.tools.tsss.momentspace.D = {fullfile(D)};
            matlabbatch{1}.spm.tools.tsss.momentspace.condthresh = 80;       

            spm_jobman('run', matlabbatch);
        end
        D = reload(D);
        
%% calculate SNR for reg par
        snr_dict = load([out_path site '-' datacat '_SNR.mat']);
        snr_idx  = strcmp(snr_dict.([site '_' datacat '_snr']).(par.meg{1,meg}).label, mfname);
        SNR      = snr_dict.([site '_' datacat '_snr']).(par.meg{1,meg}).snr(snr_idx);

        regs = eval(reg_form);
        
        %regs=5; %change it

%% Model head (Sphere model for phantom)
        clear matlabbatch
        matlabbatch{1}.spm.meeg.source.headmodel.D = {fullfile(D)};
        matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
        matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
        matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.mri = {mrifile};
        matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'Nasion';
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = [0 79.5 0];
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'LPA';
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = [-79.5 0 0];
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'RPA';
        matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = [79.5 0 0];
        matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Sphere';
        spm_jobman('run', matlabbatch);

        D = reload(D);
        D.inv{1}.forward.vol.o = [0 0 0];
        D.inv{1}.forward.vol.r = 0.075;%0.0795;
        spm_eeg_inv_checkforward(D);
        save(D);
%% Apply beamforming
    for reg=regs
        out_tmp_dir=[out_dir 'BF_' mfname '-chan_' par.meg{meg} '-sss' num2str(sss) '-runtry' num2str(par.runtry)];
        mkdir(out_tmp_dir)
        cd(out_tmp_dir)
       
        if isequal(apply_lcmv, 'yes') % Apply LCMV
            disp('Running LCMV ----------------------->')
            if reg == regs(1)
                clear matlabbatch;
                EN_SPM_run_LCMV(D, par.gridres, par.bpfreq, reg, woi, par.meg, meg)
            else
                load([out_tmp_dir '/BF_LCMV.mat'])
                unix(['rm -f ' [out_tmp_dir '/BF_LCMV.mat']])
                save([out_tmp_dir '/BF.mat'], 'data', 'sources')
                clear data source features inverse output  
                clear matlabbatch
                EN_SPM_run_LCMV2(out_tmp_dir, par.bpfreq, reg, woi, par.meg, meg)
            end
            
            % Write hotspot (dipole) location in a text file
            movefile('BF.mat', 'BF_LCMV.mat')
            BF = load('BF_LCMV.mat');
            %[hval, hind]=max(BF.output.image.val);
            [~, hind]=max(abs(BF.output.image.val));
            hval = BF.output.image.val(hind);
            hspot=BF.sources.pos(hind, :)*1000;

            difff=sqrt(sum((actual_diploc(dip,:)-hspot).^2));
            disp([hspot, hval, difff])

            resultfile=[out_dir 'SPM_BF' '-' site '_estimated-source_loc.csv'];
            if meg==chancat(1) && pp==maxf(1) && sub==amps(1) && dip==dipoles(1) && reg==regs(1)
                fid = fopen(strcat(resultfile), 'a+');
                fprintf(fid, '\nGrid res = %s,\n', num2str(par.gridres));
                fprintf(fid, 'Band pass= %s,\n', num2str(par.bpfreq));
                fprintf(fid, 'reg_form = %s,\n', reg_form);
                fprintf(fid, 'Date & time = %s,\n', datestr(now));
                fprintf(fid, '\n');
                fclose(fid);  
            else
            end
            
            fid = fopen(strcat(resultfile), 'a+');
            fprintf(fid, '%s,', mfname);
            fprintf(fid, '%s,', char(par.subjects(sub)));            
            fprintf(fid, '%f,', [dip, hspot(1),hspot(2),hspot(3), hval, difff, D.ntrials,  D.nchannels, SNR, reg]);
            fprintf(fid, '%s,', par.prep{1,pp});
            fprintf(fid, '%s,', 'LCMV');
            fprintf(fid, '\n');
            fclose(fid);            
        else
        end
        
        if isequal(apply_dics, 'yes')  % Apply DICS
            clear matlabbatch
            disp('Running DICS ----------------------->')
            EN_SPM_run_DICS(D, par.gridres, reg, woi, dics_freq, par.meg, meg)
            
            % Write hotspot (dipole) location in a text file
            movefile('BF.mat', 'BF_DICS.mat')
            BF = load('BF_DICS.mat');
            %[hval, hind]=max(BF.output.image.val);
            [~, hind]=max(abs(BF.output.image.val));
            hval = BF.output.image.val(hind);
            hspot=BF.sources.pos(hind, :)*1000;
            difff=sqrt(sum((actual_diploc(dip,:)-hspot).^2));
            disp([hspot, hval, difff]);
                        
            resultfile_dics=[out_dir 'SPM_BF' '-' site '_estimated-source_loc_DICS.csv'];
            if meg==chancat(1) && pp==maxf(1) && sub==amps(1) && dip==dipoles(1) && reg==regs(1)
                fid = fopen(strcat(resultfile_dics), 'a+');
                fprintf(fid, '\nGrid res = %s,\n', num2str(par.gridres));
                fprintf(fid, 'Band pass= %s,\n', num2str(par.bpfreq));
                fprintf(fid, 'reg_form = %s,\n', reg_form);
                fprintf(fid, 'Date & time = %s,\n', datestr(now));
                fprintf(fid, '\n');
                fclose(fid);  
            else
            end
            fid = fopen(strcat(resultfile_dics), 'a+');
            fprintf(fid, '%s,', mfname);
            fprintf(fid, '%s,', char(par.subjects(sub)));            
            fprintf(fid, '%f,', [dip, hspot(1),hspot(2),hspot(3), hval, difff, D.ntrials,  D.nchannels, SNR, reg]);
            fprintf(fid, '%s,', par.prep{1,pp});
            fprintf(fid, '%s,', 'DICS');
            fprintf(fid, '\n');
            fclose(fid);
        else
        end
        
        if reg==regs(end)           
            cd('~//spm_temp//')
            unix(['rm -rf ' out_tmp_dir])
        else
            cd(out_tmp_dir)
            unix('find . -type f -name "*.nii" -delete')
        end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close all
        clearex meg pp sss sub dip allsite chancat datacat SSS apply_dics apply_lcmv out_tmp_dir...
            reg_compare dics_freq site maxf amps dipoles phantom_moving iii actual_diploc homedir...
            reg regs D out_dir trialwin woi matlabbatch sss mfname par reg_form SNR out_path 
    toc;
    end
    cd('~//spm_temp//')
    unix(['rm -rf ' out_tmp_dir])
                    end
                end
            end
            unix('find . -type f -name  "*spmeeg_*"  -delete')
            unix('find . -type f -name  "*spmeeg_*"  -delete')
        end
    end
end
cd ~
toc;
%***********************************END************************************

