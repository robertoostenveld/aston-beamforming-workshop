function [trl]= EN_make_trialdef_for_SPMbatch5(dataset, stimchan, stimval, trialwin, minGap, rmtrialstart, rmtrialend)
    % The function is intended to make trial definition file for epoching phantom data in SPM batch
    % By amit @ Elekta Neuromag (Helsinki)
    % Usage: [trl, conditionlabels]= EN_make_trialdef_for_SPMbatch('trigchan', trigval, [timewin], mingap, 'condlabel')
    % mingap: it's time in second to avoid offset to detected as additional triggers. 
    %         It should be in between trigger length and ISI.
    % Example: [trl, conditionlabels]= EN_make_trialdef_for_SPMbatch('STI201', 7, [-0.1 0.1], 0.15, 'dip7');
    %           here trigger length is 100ms and ISI is 350ms.
    cfg                     = [];                   
    cfg.dataset             = dataset;
    cfg.trialdef.eventtype  = stimchan;
    cfg.trialdef.eventvalue = stimval;                     
    cfg.trialdef.prestim    = abs(trialwin(1));                     
    cfg.trialdef.poststim   = trialwin(2); 
    cfg.trialfun            = 'EN_phantom_trialfun';
    cfg.trig_min_gap        = minGap;
    cfg.minlength = cfg.trialdef.prestim + cfg.trialdef.poststim;
    cfg = ft_definetrial(cfg);
    % Leave the first and last trigger to avoid data deficiency error
    cfg.trl(1:rmtrialstart,:)  =[];
    cfg.trl(end-rmtrialend:end,:)=[];

    trl=[cfg.trl(:,1:2), cfg.trl(:,3)];

    fprintf(['**Subtracted first ' num2str(rmtrialstart) ' and last ' num2str(rmtrialend) ' trials to remove the reset artifacts. Resulting ' num2str(size(trl,1)) ' trials in total...........\n\n'])

end
