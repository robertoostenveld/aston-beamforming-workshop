function [actual_diploc, channelslist]=EN_load_plantom_actual_dip_chanlist(site)
    if isunix
        if isequal(site, 'Aalto') % For Vectorview
            actual_diploc=load('~/Dropbox/MATLAB/aalto_phantom.mat'); 
            actual_diploc=actual_diploc.aalto_phantom;
        elseif any(strcmp({'Aston', 'Aston_mov', 'Bari', 'Biomag'}, site)) % For Triux
            actual_diploc=load('~/Dropbox/MATLAB/biomag_phantom.mat');
            actual_diploc=actual_diploc.biomag_phantom;
        end
        channelslist=load('~/Dropbox/MATLAB/channelslist.mat');
        channelslist=channelslist.channelslist;
    elseif ispc
        if isequal(site, 'Aalto') % For Vectorview
            actual_diploc=load('C:\Users\fijaiami//Dropbox//MATLAB//aalto_phantom.mat'); 
            actual_diploc=actual_diploc.aalto_phantom;
        elseif any(strcmp({'Aston', 'Aston_mov', 'Bari', 'Biomag'}, site)) % For Triux
            actual_diploc=load('C:\Users\fijaiami//Dropbox//MATLAB//biomag_phantom.mat');
            actual_diploc=actual_diploc.biomag_phantom;
        end
        channelslist=load('C:\Users\fijaiami//Dropbox//MATLAB//channelslist.mat');
        channelslist=channelslist.channelslist;
    end
end
    
    