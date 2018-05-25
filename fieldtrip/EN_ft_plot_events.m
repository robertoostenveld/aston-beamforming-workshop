function ft_plot_events(figure,cfg, keyset, valueset)
%%
% Usage: FT_PLOT_EVENTS plots timelocked events as assigned in keyset and valueset
%        && also shows the number of trials for each event category
% Inputs: 
%           figure   = figure 
%           cfg      = cfg    , defined by ft_definetrial
%           keyset   = 1xn cell array of event names,  where n= no. of event categories
%           valueset = 1xn cell array of event values, where n= no. of event categories
% Example: 
%           figure   = figure;
%           cfg      = cfg;
%           keyset   = {'VEF-UR', 'VEF-LR', 'AEF-Le', 'VEF-LL', 'AEF-Re', 'VEF-UL', 'SEF-Lh', 'SEF-Rh'};
%           valueset = [1,2, 3, 4, 5, 8, 16, 32];
%
%           ft_plot_events(figure,cfg, keyset, valueset)
%
% Author: Amit Jaiswal @ MEGIN (Elekta Oy), Helsinki, Finland

% Create axes
    axes1 = axes('Parent',figure,'YTick',valueset);

    box(axes1,'on');
    hold(axes1,'all');

% Create plot
    plot(cfg.trl(:,1),cfg.trl(:,4),'Marker','o','LineStyle','none','Color',[1 0 0],...
        'DisplayName','stims');
    xlim([cfg.trl(1,1) cfg.trl(end,2)])
    ylabel('Event id ---->', 'fontsize', 12);
    xlabel('Latency ---->', 'fontsize', 12);
    title('Events plot', 'fontsize', 14);

    clear wrightlegs
for ii=1:length(keyset)
    wrightlegs{1, ii}=char(strcat(keyset(ii), ' = ', num2str(length(find(cfg.trl(:,4) == valueset(ii))))));
end

% Create textbox
annotation(figure,'textbox',...
    [0.91 0.71 0.1 0.1],...
    'LineWidth', 0.00001, 'String',wrightlegs,'fontsize', 12,...
    'FitBoxToText','on');

