function EN_SPM_run_LCMV2(out_tmp_dir, bpfreq, reg, woi, chan_categories, meg)
%   matlabbatch{1,3}=[];
%   matlabbatch{1,4}=[];
%   matlabbatch{1,5}=[];
%   matlabbatch{1,6}=[];
%   matlabbatch{1,7}=[];
%     clear matlabbatch
    matlabbatch{1}.spm.tools.beamforming.features.BF = {[out_tmp_dir '/BF.mat']};
    matlabbatch{1}.spm.tools.beamforming.features.whatconditions.all = 1;
    matlabbatch{1}.spm.tools.beamforming.features.woi = woi;
    if strcmp(chan_categories{meg,1}, 'all')
        matlabbatch{1}.spm.tools.beamforming.features.modality = {
                                                                  'MEG'
                                                                  'MEGPLANAR'
                                                                  }';
        matlabbatch{1}.spm.tools.beamforming.features.fuse = 'meg'; %all
    elseif strcmp(chan_categories{meg,1}, 'mag')
        matlabbatch{1}.spm.tools.beamforming.features.modality = {'MEGMAG'}; %'MEG' check this if not works
        matlabbatch{1}.spm.tools.beamforming.features.fuse = 'no';
    elseif strcmp(chan_categories{meg,1}, 'grad')
        matlabbatch{1}.spm.tools.beamforming.features.modality = {'MEGPLANAR'};
        matlabbatch{1}.spm.tools.beamforming.features.fuse = 'no';
    end
    matlabbatch{1}.spm.tools.beamforming.features.plugin.cov.foi = bpfreq;
    matlabbatch{1}.spm.tools.beamforming.features.plugin.cov.taper = 'hanning';
    matlabbatch{1}.spm.tools.beamforming.features.regularisation.manual.lambda = reg; %par.reg;
    matlabbatch{1}.spm.tools.beamforming.features.bootstrap = false;
    %%matlabbatch{2}.spm.tools.beamforming.inverse.BF(1) = cfg_dep('Covariance features: BF.mat file', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));       
    matlabbatch{2}.spm.tools.beamforming.inverse.BF(1) = cfg_dep('Covariance features: BF.mat file', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
    matlabbatch{2}.spm.tools.beamforming.inverse.plugin.lcmv.orient = true;
    matlabbatch{2}.spm.tools.beamforming.inverse.plugin.lcmv.keeplf = false;
    %%matlabbatch{3}.spm.tools.beamforming.output.BF(1) = cfg_dep('Inverse solution: BF.mat file', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
    matlabbatch{3}.spm.tools.beamforming.output.BF(1) = cfg_dep('Inverse solution: BF.mat file', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
    matlabbatch{3}.spm.tools.beamforming.output.plugin.image_power.powermethod = 'trace';
    matlabbatch{3}.spm.tools.beamforming.output.plugin.image_power.whatconditions.all = 1;
    matlabbatch{3}.spm.tools.beamforming.output.plugin.image_power.sametrials = false;
    matlabbatch{3}.spm.tools.beamforming.output.plugin.image_power.woi = woi;
    matlabbatch{3}.spm.tools.beamforming.output.plugin.image_power.foi = bpfreq;
    matlabbatch{3}.spm.tools.beamforming.output.plugin.image_power.contrast = [-1 1];
    matlabbatch{3}.spm.tools.beamforming.output.plugin.image_power.logpower = true;
    matlabbatch{3}.spm.tools.beamforming.output.plugin.image_power.result = 'singleimage'; %  bycondition
    matlabbatch{3}.spm.tools.beamforming.output.plugin.image_power.scale = 2;
    if strcmp(chan_categories{1,meg}, 'all')
        matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.modality = 'MEG';
    elseif strcmp(chan_categories{1,meg}, 'mag')
        matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.modality = 'MEG'; 
    elseif strcmp(chan_categories{1,meg}, 'grad')
        matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.modality = 'MEGPLANAR';
    end
        %%matlabbatch{4}.spm.tools.beamforming.write.BF(1) = cfg_dep('Output: BF.mat file', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
    % % matlabbatch{4}.spm.tools.beamforming.write.BF(1) = cfg_dep('Output: BF.mat file', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
    % % matlabbatch{4}.spm.tools.beamforming.write.plugin.nifti.normalise = 'separate';
    % % matlabbatch{4}.spm.tools.beamforming.write.plugin.nifti.space = 'native';%'mni';
    % % 
    % % matlabbatch{5}.spm.util.disp.data(1) = cfg_dep('Write: Output files', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));

    spm_jobman('run', matlabbatch);
    
    % saveas(gcf,['SPM_' mfname '_sss' num2str(sss) '-' num2str(par.runtry) '-LCMV'],'png');
    % overlay=ft_read_mri(['uv_pow_effspmeeg_' mfname '.nii'])
            
end