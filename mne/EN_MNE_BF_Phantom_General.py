#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 16:00:50 2017

@author: amit && Alex

==================================================================================
Compute LCMV/DICS beamformer (Volumetric) for all phantom data (moving or static)
==================================================================================
>> Beamformer comparison works

"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
import os
import mne
from nilearn.plotting import plot_stat_map
from nilearn.image import index_img
from datetime import datetime
plt.close('all')
print(__doc__)

mne.set_log_level('WARNING')
#%% Set data category, site etc
for site in ['Aalto', 'Aston', 'Aston_mov', 'Bari', 'Biomag']: 
    if site=='Aalto':
        maxf, amps, dipoles = ['', '_sss'], ['20', '100','200', '1000'], (5,6,7,8)
    if site== 'Aston':
        maxf, amps, dipoles = ['','_sss'], ['20','200','1000'], (5,6,7,8,9,10,11,12)
    if site== 'Bari':
        maxf, amps, dipoles = ['', '_sss'], ['25'], np.arange(5,13)
    if site== 'Aston_mov':
        maxf, amps, dipoles, phantom_moving = ['','_tsss_mc'], ['200'], (5,6,7,8,9,10,11,12), 'yes'
    if site== 'Biomag':
        maxf, amps, dipoles = ['', '_sss'], ['500'], np.arange(5,13)
    
    datacat = 'phantom'
    apply_lcmv='yes'
    apply_dics='yes'
    reg_compare='no'
    dics_freq=[5.0, 35.0]
    chancat=[True]#, 'grad', 'mag']
    print(maxf, amps, dipoles)
    reg_form='snr_dict[dfname][0]**4/500'

#%% Run the loop
    for meg in chancat:
        for prepcat in maxf:
            for subcat in amps:
                for dipnum in dipoles:#
                    #% Set parameters dictionary
                    par={'savefig'    : 'no',
                        'savefig_res' : '',
                        'maxfilter'   : '',
                        'rawfilter'   : '',
                        'epochsfilter': 'yes',
                        'bandpass'    : 'yes',
                        'bpfreq'      : [2, 95],
                        'calc_src_v'  : '',
                        'gridres'     : 5.0,
                        'mindist'     : 5.0,
                        'exclude'     : 0.0,
                        'bem_sol'     : 'sphere',
                        'check_trial' : '',
                        'browse'      : '',
                        'more_plots'  : '',
                        'models_plot' : '',                
                        'numtry'      : 2}
                        
                    par['event_dict']= dict(Dipole=dipnum)
                    par['event_id']= dipnum
                                
                    # Load actual dipole location for phantom
                    if site in {'Bari', 'Aston', 'Aston_mov', 'JYU', 'Biomag', 'Birmingham'}:
                        act_dip=loadmat('/net/bonsai/home/amit/Documents/MATLAB/biomag_phantom.mat') # TRIUX
                        act_dip=act_dip['biomag_phantom']
                    elif site in {'Aalto', 'Biomag_old'}:
                        act_dip=loadmat('/net/bonsai/home/amit/Documents/MATLAB/aalto_phantom.mat')  # Vectorview
                        act_dip=act_dip['aalto_phantom']
                    
                    # Phantom CT information
                    subjects_dir, subject = '/opt/freesurfer/subjects/', 'Aalto_phantom'
                    trans = subjects_dir + subject + '/mri/transforms/' + 'Aalto_phantom-trans2.fif'
                    mrifile= subjects_dir + subject + '/mri/T1.mgz'
                    surffile= subjects_dir + subject + '/bem/watershed/' + subject + '_outer_skin_surface'
                    runfile('/net/bonsai/home/amit/Dropbox/Python/EN_set_directory.py')
        #%% Set MRI/MEG data and data specific constarints
                    [data_path, fname, out_dir] = EN_set_directory(site, datacat, subcat, dipnum, prepcat)
                    dfname=os.path.split(os.path.splitext(fname)[0])[1]
                    out_path = out_dir + 'MNE//'
        
                    if site=='Aalto': # Aalto old data shared very first time 
                        raw=mne.io.read_raw_fif(fname, allow_maxshield=False, preload=True, verbose=True)
                        badch = raw.info['bads'] = []
                        if prepcat=='':
                            badch = raw.info['bads'] = ['MEG2233', 'MEG2231', 'MEG0111', 'MEG2422', 'MEG1842', 'MEG0511']  
                        par['stimch']  = 'STI201'
                        par['trialwin'] = [-0.1, 0.1]
                        par['ctrlwin']  = [-0.1, -0.0]
                        par['actiwin']  = [0.0, 0.1]
                        reject=dict(grad=4000e-13, mag=4e-12)
                        badtrial, badreason = [], ''
                        dfname=os.path.split(os.path.splitext(fname)[0])[1]
                        events = mne.find_events(raw, stim_channel=par['stimch']) 

                    if site=='Aston': # Aston shared static & moving phantom data  
                        mov=''
                        raw=mne.io.read_raw_fif(fname, allow_maxshield=False, preload=True, verbose=True)
                        badch = raw.info['bads'] = []
                        if prepcat=='':
                            badch = raw.info['bads'] = ['MEG1133', 'MEG1323', 'MEG0613', 'MEG1032'] # added 'MEG0613', 'MEG1032' 
                        par['stimch']  = 'SYS201'
                        par['trialwin'] = [-0.5, 0.5]
                        par['ctrlwin']  = [-0.5, -0.0]
                        par['actiwin']  = [0.0, 0.5]
                        reject=dict(grad=4000e-13, mag=4e-12)
                        badtrial, badreason = [0], 'resetsignal'
                        dfname=os.path.split(os.path.splitext(fname)[0])[1]
                        events = mne.find_events(raw, stim_channel=par['stimch']) # Read even
                        
                    if site=='Biomag': # Aston shared static & moving phantom data  
                        raw=mne.io.read_raw_fif(fname, allow_maxshield=False, preload=True, verbose=True)
                        badch = raw.info['bads'] = []
                        if prepcat=='':
                            badch = raw.info['bads'] = ['MEG0342',  'MEG0542', 'MEG1013'] # added 'MEG0613', 'MEG1032' 
                        par['stimch']  = 'SYS201'
                        par['trialwin'] = [-0.1, 0.1]
                        par['ctrlwin']  = [-0.1, -0.0]
                        par['actiwin']  = [0.0, 0.1]
                        reject=dict(grad=4000e-13, mag=4e-12)
                        badtrial, badreason = [0], 'resetsignal'
                        dfname=os.path.split(os.path.splitext(fname)[0])[1]
                        events = mne.find_events(raw, stim_channel=par['stimch']) # Read even
                        
                    if site=='Bari': # Bari static phantom data  
                        raw=mne.io.read_raw_fif(fname, allow_maxshield=True, preload=True, verbose=True)
                        par['stimch']  = 'SYS201'
                        dfname=os.path.split(os.path.splitext(fname)[0])[1] #+ '-dip'+str(dipnum)
                        if prepcat=='':
                            badch = raw.info['bads'] = ['MEG0943', 'MEG0222', 'MEG1522', 'MEG1512', 'MEG1432', 'MEG1113', 'MEG0631']
                        par['trialwin'] = [-0.1, 0.1]
                        par['ctrlwin']  = [-0.1, -0.0]
                        par['actiwin']  = [0.0, 0.1]
                        reject=dict(grad=4000e-13, mag=4e-12)
                        badtrial, badreason = [], 'resetsignal'
                        events = mne.find_events(raw, stim_channel=par['stimch']) # Read event
                        
                    if site=='Aston_mov': # Aston shared static & moving phantom data  
                        mov='_movement'
                        raw=mne.io.read_raw_fif(fname, allow_maxshield=False, preload=True, verbose=True)
                        badch = raw.info['bads'] = []
                        if prepcat=='':
                            badch = raw.info['bads'] = ['MEG1133', 'MEG1323', 'MEG0613', 'MEG1032'] # added 'MEG0613', 'MEG1032' 
                        par['stimch']  = 'SYS201'
                        par['trialwin'] = [-0.5, 0.5]
                        par['ctrlwin']  = [-0.5, -0.0]
                        par['actiwin']  = [0.0, 0.5]
                        reject=dict(grad=4000e-13, mag=4e-12)
                        badtrial, badreason = [0], 'resetsignal'
                        dfname=os.path.split(os.path.splitext(fname)[0])[1]
                        events = mne.find_events(raw, stim_channel=par['stimch']) # Read even
                        events[:,2]=events[:,2]-events[:,1]
                        events[:,1]=events[:,1]-events[:,1]
                     
                    snr_dict=np.load(out_dir + site + '_' + datacat + '_SNR_' + str(meg) + '.npy').item()
                    raw.pick_types(meg=meg, exclude='bads')
                    raw.info
    
                    if not os.path.exists(out_path):
                            os.mkdir(out_path)
                    resultfile_lcmv = '/net/bonsai/home/amit/Dropbox/BeamComp_Resultfiles/' + site + '_phantom_data/' 'MNE/'+ 'MNEp_Result-numtry' + str(par['numtry'])+ '-phantom-'+ site + '_source_loc.csv'
                    resultfile_dics = resultfile_lcmv[:-4]+'-DICS.csv'
                    if meg==chancat[0] and subcat==amps[0] and prepcat==maxf[0] and dipnum==dipoles[0]:
                        fid = open(resultfile_lcmv, 'a+')
                        fid.writelines('=======================================================================================\n'+
                                       'raw.info["bads"] = %s'%raw.info['bads'] + '(only for raw data not for MaxFiltered)\n' +
                                       'par["trialwin"]= %s'%par['trialwin']+ '\n' +
                                       'par["ctrlwin"]= %s'%par['ctrlwin']+ '\n' +
                                       'par["actiwin"]= %s'%par['actiwin']+ '\n' +
                                       'par["bpfreq"]= %s'%par['bpfreq']+ '\n' +
                                       'par["gridres"]= %s'%par['gridres']+ '\n' +
                                       'par["mindist"]= %s'%par['mindist']+ '\n' +
                                       'par["exclude"]= %s'%par['exclude']+ '\n' +
                                       'par["bem_sol"]= %s'%par['bem_sol']+ '\n' +
                                       'Date & Time = %s' %str(datetime.now())+ '\n' +
                                       'Regularization method= %s'%reg_form +'\n' +
                                       '=======================================================================================\n')
                        fid.close()
                        if apply_dics=='yes':
                            fid = open(resultfile_dics, 'a+')
                            fid.writelines('=======================================================================================\n'+
                                           'raw.info["bads"] = %s'%raw.info['bads'] + '(only for raw data not for MaxFiltered)\n' +
                                           'par["trialwin"]= %s'%par['trialwin']+ '\n' +
                                           'par["ctrlwin"]= %s'%par['ctrlwin']+ '\n' +
                                           'par["actiwin"]= %s'%par['actiwin']+ '\n' +
                                           'par["bpfreq"]= %s'%par['bpfreq']+ '\n' +
                                           'par["gridres"]= %s'%par['gridres']+ '\n' +
                                           'par["mindist"]= %s'%par['mindist']+ '\n' +
                                           'par["exclude"]= %s'%par['exclude']+ '\n' +
                                           'par["bem_sol"]= %s'%par['bem_sol']+ '\n' +
                                           'Date & Time = %s' %str(datetime.now())+ '\n' +
                                           'Regularization method= %s'%reg_form +'\n' +
                                           '=======================================================================================\n')
                            fid.close()
                        
        #%% Browse data
                    if par['browse']=='yes':
                        mne.viz.plot_events(events, sfreq=None, first_samp=0, color=None, event_id=par['event_dict'], 
                            axes=None, equal_spacing=True, show=True)
                        raw.plot(events=events, duration=3.0, start=0.0, n_channels=20, bgcolor='w', color=None,
                                 bad_color=(0.5, 0.1, 0.1), event_color='cyan', scalings=None, remove_dc=True,
                                 order='type', show_options=False, title=fname, show=True, block=False,
                                 highpass=None, lowpass=None, filtorder=4, clipping=None, show_first_samp=False)
                        raw.plot_psd(tmin=0.0, tmax=raw.times[-1], fmin=0, fmax=300, proj=False, n_fft=None, picks=None, ax=None, 
                                    color='black', area_mode='std', area_alpha=0.33, n_overlap=0, dB=True, average=False,
                                    show=True, n_jobs=1, line_alpha=0.5, spatial_colors=True, xscale='linear', verbose=True)
        
        #%% Apply filter if required
                    if par['maxfilter']=='yes'and prepcat=='':
                        raw.fix_mag_coil_types()
                        raw=mne.preprocessing.maxwell_filter(raw, origin=(0., 0., 0.0), int_order=8, ext_order=3, 
                                calibration=None, cross_talk=None, st_duration=None, st_correlation=0.98, coord_frame='head',
                                destination=None, regularize='in', ignore_ref=False, bad_condition='error', head_pos=None,
                                st_fixed=True, st_only=False, mag_scale=100.0, verbose=None)
                    if par['rawfilter']=='yes' and par['epochsfilter']!='yes':
                        raw.notch_filter(50, filter_length='auto', phase='zero') # Notch filter
                        raw.filter(par['bpfreq'][0], par['bpfreq'][1], l_trans_bandwidth=min(max(2 * 0.01, 2), 2), 
                                   h_trans_bandwidth=min(max(70 * 0.01, 2.), raw.info['sfreq'] / 2. - 70), 
                                   filter_length='auto', phase='zero')

        #%% Epoch data
                    raw.info.normalize_proj()
                    epochs = mne.Epochs(raw, events, par['event_id'], par['trialwin'][0], 
                                        par['trialwin'][1], baseline=tuple(par['ctrlwin']), picks=None, 
                                        preload=True, reject=None, flat=None, proj=False, decim=1,
                                        reject_tmin=None, reject_tmax=None, detrend= None, 
                                        on_missing='error', reject_by_annotation=True, verbose=True)
                    epochs.drop(badtrial, badreason)
                    if site=='Aston' and prepcat=='_sss' and subcat=='200' and dipnum==12:
                        epochs.drop([37,38,39]) 
                    del raw
                    
                    if par['epochsfilter']=='yes' and par['rawfilter']!='yes':
                        epochs.filter(par['bpfreq'][0], par['bpfreq'][1], picks=None, filter_length=100, l_trans_bandwidth=2, 
                                      h_trans_bandwidth=2, n_jobs=1, method='fir', iir_params=None, 
                                      phase='zero', fir_window='hamming', fir_design='firwin2', pad='edge', 
                                      verbose=True)
                        epochs.filter(51, 49, picks=None, filter_length=100, l_trans_bandwidth=2, 
                                      h_trans_bandwidth=2, n_jobs=1, method='fir', iir_params=None, 
                                      phase='zero', fir_window='hamming', fir_design='firwin2', pad='edge', 
                                      verbose=True)
                    
                    if par['check_trial']=='yes':
                        epochs.plot(picks=None, scalings=None, n_epochs=10, n_channels=30, event_colors=None, 
                                    title='Epochs plot (cascaded)', events=None, show=True, block=False)
                        
        #%% Average & plot            
                    evoked= epochs.average()

        #%% Model BEM
                    bem=mne.make_sphere_model(r0=(0.0, 0.0, 0.0), head_radius=None, info=None, verbose=True)
        
        #%% Compute Source Space .................................>>>>
                    if par['calc_src_v']=='yes':
                        src_v=mne.setup_volume_source_space(subject=None, pos=par['gridres'], mri=mrifile,
                                    bem=None, surface=surffile, mindist=par['mindist'], exclude=par['exclude'], 
                                    subjects_dir=None, volume_label=None, add_interpolator=True, verbose=True)
                        print(src_v)
                        mne.write_source_spaces(subjects_dir + subject + '/' + subject + '-volumetric-gridres' + str(par['gridres'])+ 'mm-mindist'+ str(par['mindist'])+ 'mm-exclude'+ str(par['exclude']) + 'mm-src.fif', src_v, overwrite=True)
                    else:
                        src_v=subjects_dir + subject + '/' + subject + '-volumetric-gridres' + str(par['gridres'])+ 'mm-mindist'+ str(par['mindist'])+ 'mm-exclude'+ str(par['exclude']) + 'mm-src.fif'
                        src_v=mne.read_source_spaces(src_v, patch_stats=False, verbose=True)
                        
                    src_v=subjects_dir + subject + '/' + subject + '-volumetric-gridres' + str(par['gridres'])+ 'mm-mindist'+ str(par['mindist'])+ 'mm-exclude'+ str(par['exclude']) + 'mm-src.fif'

        #%% Forward solution..................>> 
                    fwd = mne.make_forward_solution(epochs.info, trans=trans, src=src_v, bem=bem, 
                                                    meg=meg, eeg=False, mindist=None, n_jobs=1)
                    leadfield = fwd['sol']['data']
                    print("Leadfield size : %d sensors x %d dipoles" % leadfield.shape)
        
        #%% Apply time domain Beamforming (LCMV):
                    if apply_lcmv=='yes':
                        reg=eval(reg_form)
                        noise_cov = mne.compute_covariance(epochs, tmin=par['ctrlwin'][0], tmax=par['ctrlwin'][1], method='empirical') 
                        data_cov = mne.compute_covariance(epochs, tmin=par['actiwin'][0], tmax=par['actiwin'][1], method='empirical')   

                        filters = mne.beamformer.make_lcmv(evoked.info, fwd, data_cov, reg=reg,noise_cov=noise_cov,
                                                                  pick_ori='max-power', 
                                                                  rank=np.sum(np.linalg.svd(noise_cov.data)[1]>1e-25),
                                                                  weight_norm='unit-noise-gain', reduce_rank=True)
                        stc = mne.beamformer.apply_lcmv(evoked, filters, max_ori_out='signed')
                        
                        stc=np.abs(stc)

                        src_peak, t_peak=stc.get_peak()
                        est_loc = fwd['src'][0]['rr'][src_peak]*1000 
                        loc_err=np.sqrt(np.sum(np.square(act_dip[dipnum-1]-est_loc))); # Calculate loc_err
        
                        print('Act_Sourceloc for %s' %dfname + '= %s' % str(act_dip[dipnum-1])) 
                        print('Est_SourceLoc for %s' %dfname + '= %s' % str(np.around(est_loc,1)))
                        print('Peak_Value for %s' %dfname + '= %.2f' % stc.data.max())
                        print('Loc_error for %s' %dfname + '= %.1f' % loc_err)

                        if par['more_plots']:
                            plt.figure()
                            ts_show = -50  
                            plt.plot(1e3 * stc.times, stc.data[np.argsort(stc.data.max(axis=1))[ts_show:]].T)
                            plt.title(dfname + ' for %d largest sources'%abs(ts_show))
                            plt.xlabel('time (ms)')
                            plt.ylabel('%s value ' %'LCMV stc'+ '@reg=%.2f'%reg )
                            plt.show()
                        
                            stc.crop(0.0, stc.times[-1])
                            thresh = stc.data.max()*2/100
                            timepoint = int(t_peak//stc.tstep - stc.times[0]//stc.tstep)
        
                            img=mne.save_stc_as_volume('lcmv_stc.nii', stc, fwd['src'], dest='mri', mri_resolution=False)
                            plot_stat_map(index_img(img,  int(timepoint)), 
                                          mrifile, draw_cross=True, threshold=thresh, 
                                          title = dfname +  '/ LCMV/ '+ 'reg=%.2f'%reg + '/ Vpeak=%.3f\n'%stc.data.max()+ 
                                          'Tpeak=%.3fs'%t_peak + '/ Est_loc= [%.1f, %.1f, %.1f]' %tuple(est_loc) + 
                                          '/ Loc_err= %.1f' % loc_err) 
                            os.remove('lcmv_stc.nii')
                                                
                        fid = open(resultfile_lcmv, 'a+')
                        fid.writelines('%s,' %dfname)
                        fid.writelines('%s,' %subcat)
                        fid.writelines('%.f,' %dipnum)
                        fid.writelines('%.2f,' %est_loc[0])
                        fid.writelines('%.2f,' %est_loc[1])
                        fid.writelines('%.2f,' %est_loc[2])
                        fid.writelines('%.2f,' %stc.data.max())
                        fid.writelines('%.2f,' %loc_err)
                        fid.writelines('%.f,' %len(evoked.ch_names))
                        fid.writelines('%.3f,' %snr_dict[dfname][0])
                        fid.writelines('%.3f,' %reg)
                        fid.writelines('%s,' %prepcat)
                        fid.writelines('%s\n' %'LCMV')
                        fid.close()
                        
                        plt.close('all')
        #%% Compare the effect of regilarization parameter 
                        if reg_compare=='yes':
                            trial=1
                            regs=[0.001, 0.005, 0.01, 0.05,0.1,0.5,1,2,4,8,10,20,50,75,100,150]
                            if par['savefig_res']=='yes':
                                plt.ion()  # make plot interactive
                                _, ax = plt.subplots(4, 4)  # create subplots
                                plt.tight_layout()
                                ts_show = -50  # show the 200 largest responses
                                plt.suptitle('LCMV STC plot for %s-'%dfname + ' for %d largest sources '%abs(ts_show)+ 'BEM(sphere model)', fontsize=18)
                            
                            cnt=0
                            for reg in regs:
                                cnt=cnt+1
                                filters = mne.beamformer.make_lcmv(evoked.info, fwd, data_cov, reg=reg,noise_cov=noise_cov,
                                                                  pick_ori='max-power', 
                                                                  rank=np.sum(np.linalg.svd(noise_cov.data)[1]>1e-25),
                                                                  weight_norm='unit-noise-gain', reduce_rank=True)
                                stc = mne.beamformer.apply_lcmv(evoked, filters, max_ori_out='signed')
                                stc=np.abs(stc)
                                v_peak, t_peak=stc.get_peak()
                                est_loc=fwd['src'][0]['rr'][v_peak]*1000 # in mili meter
                                loc_err=np.sqrt(np.sum(np.square(act_dip[dipnum-1]-est_loc))); # Calculate loc_err
                                
                                regcompresfile=resultfile_lcmv[:-4]+'-Reg_Compare.csv'
                                fid = open(regcompresfile, 'a+')
                                fid.writelines('%s,' %dfname)
                                fid.writelines('%s,' %subcat)
                                fid.writelines('%.f,' %dipnum)
                                fid.writelines('%.2f,' %est_loc[0])
                                fid.writelines('%.2f,' %est_loc[1])
                                fid.writelines('%.2f,' %est_loc[2])
                                fid.writelines('%.2f,' %stc.data.max())
                                fid.writelines('%.2f,' %loc_err)
                                fid.writelines('%.f,' %len(evoked.ch_names))
                                fid.writelines('%.f,' %evoked.nave)
                                fid.writelines('%.2f,' %reg)
                                fid.writelines('%s,' %prepcat)
                                fid.writelines('%s\n' %'LCMV')
                                fid.close()
                        
                                if par['savefig_res']=='yes':
                                    ax = plt.subplot(4,4,cnt)
                                    plt.plot(1e3 * stc.times, stc.data[np.argsort(stc.data.max(axis=1))[ts_show:]].T)
                                    #ax.set_xlabel('time (ms)')
                                    ax.set_title('EstLoc=%s ' %np.round(est_loc,2) + '/PeakVal=%.1f ' %stc.data.max() +
                                                 '/LocErr=%.1f ' % loc_err + '/@%dms'%(t_peak*1000), fontsize=11, y=0.99)
                                    ax.set_ylabel('%s ' % 'LCMV' + '@Reg= %.3f' %reg)
                                    print(str(reg) + '  done >>>>>>>>>>>>>>>')
                                    print('Loc_error for %s' %dfname + '= %.1f' % loc_err)
                                    
                            plt.pause(0.25)
                            manager=plt.get_current_fig_manager()
                            manager.window.showMaximized()
                            plt.pause(0.25)
                            plt.tight_layout()
                            plt.pause(0.25)
                            plt.subplots_adjust(top=0.94, bottom=0.01)
                            plt.pause(0.25)
                            plt.subplots_adjust(hspace=0.22, wspace=0.16)
                            figname=out_path + 'Regularization_effects_LCMV_' + dfname + str(trial) + '.png'
                            while os.path.exists(figname):
                                trial=trial+1
                                figname=out_path + 'Regularization_effects_LCMV_' + dfname + str(trial) + '.png'
                            plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, 
                                                format=None, transparent=False, bbox_inches='tight', pad_inches=0.2, frameon=None)
                            plt.close()
        #%% Apply frequency domain beamforming (DICS):
                    if apply_dics=='yes': 
                        noise_csd = mne.time_frequency.csd_epochs(epochs, mode='multitaper', tmin=par['ctrlwin'][0], 
                                                                  tmax=par['ctrlwin'][1], fmin=dics_freq[0], fmax=dics_freq[1])
                        data_csd = mne.time_frequency.csd_epochs(epochs, mode='multitaper', tmin=par['actiwin'][0], 
                                                                 tmax=par['actiwin'][1], fmin=dics_freq[0], fmax=dics_freq[1])
                
                        stc = mne.beamformer.dics(evoked, fwd, noise_csd, data_csd, reg=reg, verbose=True)
                        
                        stc=np.abs(stc)
                        src_peak, t_peak=stc.get_peak()
                        est_loc = fwd['src'][0]['rr'][src_peak]*1000 
                        loc_err=np.sqrt(np.sum(np.square(act_dip[dipnum-1]-est_loc))); # Calculate loc_err
        
                        print('Act_Sourceloc for %s' %dfname + '= %s' % str(act_dip[dipnum-1])) 
                        print('Est_SourceLoc for %s' %dfname + '= %s' % str(np.around(est_loc,1)))
                        print('Peak_Value for %s' %dfname + '= %.2f' % stc.data.max())
                        print('Loc_error for %s' %dfname + '= %.1f' % loc_err)
                        
                        if par['more_plots']:
                            plt.figure()
                            ts_show = -50  
                            plt.plot(1e3 * stc.times, stc.data[np.argsort(stc.data.max(axis=1))[ts_show:]].T)
                            plt.title(dfname + ' for %d largest sources'%abs(ts_show))
                            plt.xlabel('time (ms)')
                            plt.ylabel('%s value ' %'DICS stc'+ '@reg=%.2f'%reg )
                            plt.show()
                        
                            stc.crop(0.0, stc.times[-1])
                            thresh = stc.data.max()*2/100
                            timepoint = int(t_peak//stc.tstep - stc.times[0]//stc.tstep)
        
                            img=mne.save_stc_as_volume('dics_stc.nii', stc, fwd['src'], dest='mri', mri_resolution=False)
                            plot_stat_map(index_img(img,  int(timepoint)), 
                                          mrifile, draw_cross=True, threshold=thresh, 
                                          title = dfname +  '/ DICS/ '+ 'reg=%.2f'%reg + '/ Vpeak=%.3f\n'%stc.data.max()+ 
                                          'Tpeak=%.3fs'%t_peak + '/ Est_loc= [%.1f, %.1f, %.1f]' %tuple(est_loc) + 
                                          '/ Loc_err= %.1f' % loc_err) 
                            os.remove('lcmv_stc.nii')
                            
                        fid = open(resultfile_dics, 'a+')
                        fid.writelines('%s,' %dfname)
                        fid.writelines('%s,' %subcat)
                        fid.writelines('%.f,' %dipnum)
                        fid.writelines('%.2f,' %est_loc[0])
                        fid.writelines('%.2f,' %est_loc[1])
                        fid.writelines('%.2f,' %est_loc[2])
                        fid.writelines('%.2f,' %stc.data.max())
                        fid.writelines('%.2f,' %loc_err)
                        fid.writelines('%.f,' %len(evoked.ch_names))
                        fid.writelines('%.3f,' %snr_dict[dfname][0])
                        fid.writelines('%.3f,' %reg)
                        fid.writelines('%s,' %prepcat)
                        fid.writelines('%s\n' %'DICS')
                        fid.close()
                        
                        plt.close('all')
#%%####################################### END (Good luck)####################################################
