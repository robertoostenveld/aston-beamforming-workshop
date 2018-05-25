function [data_path, fname, out_path]=EN_set_directory2(site, datacat, subjects, sub, dip, prep, pp)
% Author: Amit Jaiswal @ MEGIN (Elekta Oy), Helsinki
% This function is to set the data directory, MEG filename  etc.

    if isunix % *************************** For Linux
        if isequal(site, 'Aalto') && isequal(datacat, 'phantom')
            data_path='/neuro/data/beamformer_data/Beamformer_share/phantom/'; 
            fname=[data_path subjects{1,sub} 'nAm/dip0' num2str(dip) '_' subjects{1,sub} 'nAm' prep{1,pp} '.fif'];
            out_path=['/net/bonsai/home/amit/Dropbox/BeamComp_Resultfiles/' site  '_phantom_data/'];

        % For Aston phantom data    
        elseif isequal(site, 'Aston') && isequal(datacat, 'phantom')
            data_path='/neuro/data/BeamComp_data_share/Aston_phantom_data/';
            fname =[data_path 'Amp' subjects{1,sub} '_IASoff/' 'Amp' subjects{1,sub} '_Dip' num2str(dip) '_IASoff' prep{1,pp} '.fif'];
            out_path=['/net/bonsai/home/amit/Dropbox/BeamComp_Resultfiles/' site  '_phantom_data/'];
            
        elseif isequal(site, 'Aston_mov') && isequal(datacat, 'phantom')
            data_path='/neuro/data/BeamComp_data_share/Aston_phantom_data/';
            fname =[data_path 'Amp' subjects{1,sub} '_IASoff_movement/' 'Amp' subjects{1,sub} '_Dip' num2str(dip) '_IASoff_movement' prep{1,pp} '.fif'];
            out_path=['/net/bonsai/home/amit/Dropbox/BeamComp_Resultfiles/' site  '_phantom_data/'];
             
        % For Bari phantom data      
        elseif isequal(site,'Bari') && isequal(datacat, 'phantom')            
%             data_path='/media/ELEKTA_HDD/DATA/Bari_visit2/Bari_data/Phantom_data/171206/'; 
%             fname=[data_path 'All_dipole_25nAmp1' prep{1,pp} '.fif'];
            data_path='/media/ELEKTA_HDD/DATA/Bari_visit2/Bari_data/Phantom_data/171206/25nAm_IASon/'; 
            fname = [data_path '25nAmp_IASon_dip' num2str(dip) prep{1,pp} '.fif'];
            out_path=['/net/bonsai/home/amit/Dropbox/BeamComp_Resultfiles/' site  '_phantom_data/']; 
        
        elseif  isequal(site,'Biomag') && isequal(datacat, 'phantom')
            data_path = '/media/ELEKTA_HDD/DATA/phantom/Biomag/161229/500nAmp_IASoff/' ;
            fname     = [data_path '500nAmp_IASoff_dip' num2str(dip) prep{1,pp} '.fif'];
            out_path  = ['/net/bonsai/home/amit/Dropbox/BeamComp_Resultfiles/' site  '_phantom_data/'];
        end

    elseif ispc %&& ~exist('D:\DATA\', 'dir') % ******************** For Windows
        main_dir='\\172.16.51.17\data\rd\ChildBrain\BeamComp//';
        main_dir2='F:\Aston_phantom_data//';
        
        % For Aalto phantom data  
        if isequal(site, 'Aalto')
            %data_path=[main_dir 'Data_master//Aalto_phantom_data//']; 
            data_path= 'C:\DATA\AV2\Beamformer_share\phantom//';
            fname=[data_path subjects{1,sub} 'nAm//dip0' num2str(dip) '_' subjects{1,sub} 'nAm' prep{1,pp} '.fif'];
            out_path=['C:\Users\fijaiami\Dropbox\BeamComp_Resultfiles//' site  '_phantom_data//'];

        % For Aston phantom data    
        elseif isequal(site, 'Aston')
            %data_path=[main_dir 'Data_master//Aston_phantom_data//'];
            data_path=main_dir2;
            fname =[data_path 'Amp' subjects{1,sub} '_IASoff/' 'Amp' subjects{1,sub} '_Dip' num2str(dip) '_IASoff' prep{1,pp} '.fif'];
            out_path=['C:\Users\fijaiami\Dropbox\BeamComp_Resultfiles//' site  '_phantom_data//'];
        
        % For Aston moving phantom data
        elseif isequal(site, 'Aston_mov')
            %data_path=[main_dir 'Data_master//Aston_phantom_data//'];
            data_path=main_dir2;
            fname =[data_path 'Amp' subjects{1,sub} '_IASoff_movement//' 'Amp' subjects{1,sub} '_Dip' num2str(dip) '_IASoff_movement' prep{1,pp} '.fif'];
            out_path=['C:\Users\fijaiami\Dropbox\BeamComp_Resultfiles//' site  '_phantom_data//'];
        
        % For Biomag phantom data  
        elseif isequal(site,'Biomag')            
            data_path='/neuro/data/phantom/biomag/161229/';
            fname=[data_path 'Alldipoles_500nAm_raw' prep{1,pp} '.fif'];
            out_path=['C:\Users\fijaiami\Dropbox\BeamComp_Resultfiles//' site  '_phantom_data//'];

        % For Bari phantom data      
        elseif isequal(site,'Bari')            
%             data_path = 'F:\';
%             fname = [data_path 'All_dipole_25nAmp1' prep{1,pp} '.fif'];
            data_path = 'F:\25nAm_IASon\';
            fname = [data_path '25nAmp_IASon_dip' num2str(dip) prep{1,pp} '.fif'];
            out_path = ['C:\Users\fijaiami\Dropbox\BeamComp_Resultfiles//' site  '_phantom_data//'];
        end
    end
end