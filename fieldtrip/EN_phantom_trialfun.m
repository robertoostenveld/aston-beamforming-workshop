function [trl, event] = EN_phantom_trialfun(cfg)
% The trial definition function for phantom data
% By amit jaiswal @ Elekta Neuromag (Helsinki)

hdr   = ft_read_header(cfg.dataset, 'checkmaxfilter', 0); %read header info
event = EN_findmytrigger(cfg.dataset, cfg.trialdef.eventtype, cfg.trialdef.eventvalue, cfg.trig_min_gap); % read events

% search for "trigger" events
value  = [event.value]';
sample = [event.sample]';

% define pretim and poststim time duration
prestim  = cfg.trialdef.prestim  * hdr.Fs;
poststim =  cfg.trialdef.poststim * hdr.Fs;

stim_indx  = find(ismember(value,cfg.trialdef.eventvalue));

trl = [];
for istim = 1:length(stim_indx);
  newtrl   = [ sample(stim_indx(istim))-prestim, sample(stim_indx(istim))+poststim, -prestim, value(stim_indx(istim))];
  trl      = [ trl; newtrl ];
end

%% Subfunction to identify triggers
function [event]=EN_findmytrigger(megfile, trigchan, trigval, triglen)
hdr=ft_read_header(megfile, 'checkmaxfilter', 0);
dat=ft_read_data(megfile, 'checkmaxfilter', 0);
ind= ismember(hdr.label,trigchan);
[x,y]=size(dat);
if y>=x
    dat=dat';
else
end
eventchan=dat(:,ind);

[a,~,~]=find(eventchan==trigval);

onset=a(1);
for ii=2:length(a)
    if a(ii)>=onset(end) + triglen * hdr.Fs
        onset(length(onset)+1)=a(ii);
    else
    end
end
onset=onset';
for jj=1:length(onset)
    event(jj,1).type=trigchan;
    event(jj,1).sample= onset(jj);
    event(jj,1).value= trigval;
    event(jj,1).offset=[];
    event(jj,1).duration=[];
end
end
end