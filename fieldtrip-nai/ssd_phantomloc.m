% Sarang's script for analysis of phantom data from Aston beamformer
% workshop (Feb 2017)
%
% dependencies: https://github.com/fieldtrip/fieldtrip [as of 27.Feb.2017]
%               https://github.com/meeg-cfin/nemolab  [contains ft_apply_spatialfilter.m for commonweights option]
%               for plotting results: SPM8/12 and https://github.com/fieldtrip/fieldtrip/contrib/nutmegtrip

polaritytweak = 0;  % = 1 to flip voxel polarities for a nicer plot; con: requires memory/time
plotresults = 0; % = 1 to pause script and plot result of each localization
badchans = {'-MEG2233', '-MEG2422','-MEG0111'};
channel =  {'MEGMAG' badchans{:}};
coilacc = 2;
commonweights = 1; % = 1 calculates weights based on 100nAm dataset and applies to all

% nemo_ftsetup % sets up paths to FieldTrip, NutmegTrip, etc.

diploc = [32.5 0 56.3
          27.5 0 47.6
          22.5 0 39.0
          17.5 0 30.3];


datasets = {'dip05_20nAm.fif','dip05_100nAm.fif','dip05_200nAm.fif','dip05_1000nAm.fif'
    'dip06_20nAm.fif','dip06_100nAm.fif','dip06_200nAm.fif','dip06_1000nAm.fif'
    'dip07_20nAm.fif','dip07_100nAm.fif','dip07_200nAm.fif','dip07_1000nAm.fif'
    'dip08_20nAm.fif','dip08_100nAm.fif','dip08_200nAm.fif','dip08_1000nAm.fif' };




%% lead field calculation
if(exist('phantomlf.mat','file'))
    load phantomlf.mat
else
    cfg = [];
    cfg.channel = channel;
    cfg.dataset=datasets{1,1}; % all datasets have same sensor positions; just read from first dataset
    cfg.coilaccuracy = coilacc;
    data = ft_preprocessing(cfg);
    
    vol = [];
    vol.o = [0 0 0];
    vol.r = 100;
    vol.unit = 'mm';
    
    cfg = [];
    cfg.grad = data.grad;
    cfg.method = 'singlesphere';
    cfg.reducerank = 'no';
    %    [x,y,z]=ndgrid(10:2:40,-10:2:10,30:2:70);
    [x,y,z]=ndgrid(10:1:40,-10:1:10,22:1:65);
    
    %cfg.grid.resolution = 5;
    cfg.grid.pos = [x(:) y(:) z(:)];
    cfg.grid.unit = 'mm',
    cfg.vol = vol;
    
    leadgrid = ft_prepare_leadfield(cfg);
    save phantomlf leadgrid vol
end


%%
for jj=1:size(datasets,1) % across dipole locations
    for ii=1:size(datasets,2) % across dipole strengths
        cfg = [];
        cfg.dataset=datasets{jj,ii}; % for weight computation
        
        % mark triggers (ft_definetrial doesn't work on this data)
        cfg.channel='STI201';
        trigdata=ft_preprocessing(cfg);
        trigs = find(diff(trigdata.trial{1})>1)';
        
        cfg.trl = [trigs-100 trigs+300];
        cfg.trl(:,3) = -100; % 1000 Hz sampling
        cfg.trl(end,:) = []; % last one is too late
        
        cfg.demean = 'yes';
        cfg.lpfilter = 'yes'; cfg.lpfreq = 40;
        
        cfg.channel = channel;
        cfg.coilaccuracy = coilacc;
        
        data = ft_preprocessing(cfg)
        data.grad = ft_convert_units(data.grad,'mm');
        
        %%
        cfgtl=[];
        cfgtl.covariance       = 'yes';
        cfgtl.covariancewindow = 'all';  % may need to change if not desirable to include pre-stim interval
        
        timelockbp{ii} = ft_timelockanalysis(cfgtl,data);
        
    end
    
    
    %%
    cfg                   = [];
    cfg.grid              = leadgrid; % leadfield, which has the grid information
    cfg.vol               = vol; % volume conduction model (headmodel) <-- FIXME: ft_sourceanalysis insists on this even if not necessary (i.e., grid already computed)
    cfg.keepfilter        = 'yes';
    cfg.method            = 'lcmv';
    cfg.(cfg.method).reducerank   = 'no';
    cfg.(cfg.method).fixedori = 'yes';
    cfg.(cfg.method).projectnoise = 'yes';
    cfg.(cfg.method).weightnorm   = 'nai'; %% NOTE: nai or lfnorm seems crucial for good performance!
    cfg.(cfg.method).keepfilter   = 'yes';
    if(commonweights) % calculate weights with 100 nAm phantom, apply to others
        source_bp{2}         = ft_sourceanalysis(cfg, timelockbp{2}); % high-frequencies only
        for ii=[1 3 4]
            source_bp{ii} = ft_apply_spatialfilter(timelockbp{ii},source_bp{2});
        end
    else % calculate weights with respective dataset (i.e., 'normal' way)
        for ii=1:length(timelockbp)
            source_bp{ii} = ft_sourceanalysis(cfg, timelockbp{ii}); % high-frequencies only
        end
    end
    %% quantify performance
    
    for ii=1:length(timelockbp)
        tmp=cell2mat(source_bp{ii}.avg.mom);
        [peakval,peakind] = max(abs(tmp(:)));
        [peakvox_idx,peaktime_idx] = ind2sub(size(tmp),peakind);
        
        diprecon{ii,jj} = source_bp{ii}.pos(peakvox_idx,:);
        dipampl(ii,jj) = peakval;
        
        diperror(ii,jj) = norm(diprecon{ii,jj} - diploc(jj,:));
    end
    
    %% plot source maps, if desired -- pauses for each image
    if(plotresults)
        for ii=1:length(timelockbp)
            source_bp{ii}.coordsys = 'spm';
            source_bp{ii}.unit = 'mm';
            
            %% adjusts polarities for nicer plotting, but takes time for dense voxel grids
            if(polaritytweak)
                cfgpol = [];
                cfgpol.toilim = [0 0.1];
                source_bp{ii} = nmt_polaritytweak(cfgpol, source_bp{ii});
            end
            
            %%
            cfgplot = [];
            cfgplot.mripath = '/Users/sarang/lab/projects/sss/phantom/CT/Elekta_Vectorview_phantom_ct.nii';
            cfgplot.funparameter = 'avg.mom';
            
            nmt_sourceplot(cfgplot,source_bp{ii});
            disp('press any key to plot next result...')
            pause
        end
    end
end
