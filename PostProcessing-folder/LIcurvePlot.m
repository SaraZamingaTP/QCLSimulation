%% Function to represent the LIV curve 
function LIcurvePlot(isExpComparison)

FileList = dir(fullfile(pwd, '*-StateVersion.mat'));
FileList2 = dir(fullfile(pwd, 'LI_*.mat'));

figLI=figure();

for FileIndex=1:length(FileList)

    %Load file
    load(FileList(FileIndex).name);

    load(FileList2(FileIndex).name);

    Rfront=str2double(extractBetween(FileList2(FileIndex).name,'Rfront','-'));
    Rback=str2double(extractBetween(FileList2(FileIndex).name,'Rback','.mat'));
    %% Generation of the plot of Power vs Current

    %LegendString=sprintf('$R_{front}$=%.2f, $R_{back}$=%.2f', Rfront, Rback);
      
    for CurrentIndex=1:length(InputCurrent)
        pl(CurrentIndex)=plot(InputCurrent(CurrentIndex), PLtimeAv(CurrentIndex), ...
        'Color', ColorState{CurrentIndex}, 'Marker', 'square', ...
        'MarkerFaceColor',ColorState{CurrentIndex});    
        %'Color', ColorState{CurrentIndex}, 'LineStyle','-', 'LineWidth', 1.5, ...
            
        hold on;
            
    end

end

%% Experimental Data

if isExpComparison
    ExpCurrent = 9:1:420; %mA
    
    ExpPowerUT = zeros(1, 283); %till 291 mA INCLUDED; %power under threshold
    
    ExpPowerAT = [0.04 0.14 0.25 0.35 0.45 0.56 0.67 0.79 0.91 ...
        1.02 1.14 1.26 1.37 1.49 1.60 1.71 1.83 1.95 2.04 2.16 2.27 2.39 2.50 2.61 2.72 2.84 2.95 3.05 3.16 ...
        3.27 3.38 3.48 3.59 3.69 3.80 3.90 4.01 4.11 4.23 4.32 4.43 4.53 4.64 4.74 4.85 4.95 5.05 5.13 5.25 ...
        5.35 5.46 5.56 5.66 5.76 5.86 5.95 6.05 6.13 6.21 6.29 6.36 6.43 6.50 6.58 6.64 6.71 6.78 6.84 6.92 ...
        7.00 7.07 7.14 7.21 7.28 7.34 7.40 7.45 7.50 7.57 7.62 7.68 7.73 7.81 7.86 7.93 7.99 8.06 8.12 8.17 ...
        8.23 8.28 8.33 8.37 8.42 8.47 8.53 8.58 8.64 8.68 8.72 8.78 8.83 8.88 8.93 8.97 9.01 9.04 9.08 9.57 ... 
        9.61 9.66 9.70 9.74 9.78 9.81 9.84 9.86 9.90 9.94 9.97 10.00 10.03 10.05 10.07 10.08 10.10 10.13 10.14 10.16 ...
        ]; %mW
    
    ExpPower = [ExpPowerUT ExpPowerAT]*100/37;

    hold on
    pl(length(FileList)+1)=plot(ExpCurrent,ExpPower, 'linewidth', 1.5, ...
        'DisplayName', 'Experiment'); 
end

%% Figure settings

% axis settings
xlabel('Current (mA)','fontsize',30,'interpreter','latex');
ylabel('Output power (mW)','fontsize',30,'interpreter','latex');

set(gca,'xminortick', 'on','yminortick', 'on','FontSize', 25,...
        'TickLabelInterpreter','latex', 'TickDir', 'in') %Ticks
xlim([290 max(InputCurrent)])
ylim([0 460])

%legend settings
%StatesLegend = {'OFF', 'SM', 'H', 'D'};
ColorLegend={'#120C5F'; '#FF7800'; '#BA1EF5'; '#1EF525'};
for ColorIndex=length(ColorLegend):-1:1
    MatchIndex(ColorIndex)=find(strcmp(ColorLegend{ColorIndex}, ColorState), 1);
end


legend(pl(MatchIndex), 'location', 'northwest','NumColumns', 1,'fontsize',25,'interpreter','latex')
%box off
%legend boxoff
      
% Size
ScreenSize=get(groot,'ScreenSize');
ScreenSize=ScreenSize(3:4);
Height=400;
Width=600;

set(gcf, 'color', 'w') ;
set(figLI,'GraphicsSmoothing','on')
figLI.Position=[(ScreenSize(1)-Width)/2 (ScreenSize(2)-Height)/2 Width Height];

end