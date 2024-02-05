function SurfOpticalSpectrumMap(SweepParameter, TStart, TEnd, isMap)

%% Optical Spectrum surf function 

%% Figure settings
ScreenSize=get(groot,'ScreenSize');
ScreenSize=ScreenSize(3:4);
Height=600;
Width=1200;

%% % Load your data into Matlab
Subsample=1;
FilenamePrefix = 'Res_'; 
% use the dir function to find all files that match the pattern
FileList = dir(fullfile(pwd, [FilenamePrefix, '*']));

%for loop to save the sweeping vector
for Index=length(FileList):-1:1
    if strcmp(SweepParameter, 'Current')==1
        SweepParam(Index,1)=str2double(extractBetween(FileList(Index).name,'_','mA'));
    elseif strcmp(SweepParameter, 'LEF')==1
        SweepParam(Index,1)=str2double(extractBetween(FileList(Index).name,'alphaf','_'));
    elseif strcmp(SweepParameter, 'kDFB')==1
        SweepParam(Index,1)=str2double(extractBetween(FileList(Index).name,'kcoup','_'));
    end
end

%sort in ascendent order the SweepParam values
SweepParam=sort(SweepParam);

%values of SweepParam I want to represent in the map
SweepParamValues=290:700;
% Create a logical index that checks if each element of array1 is in array2
LogicalIndex = ismember(SweepParam, SweepParamValues);

% Use the logical index to extract the elements from array1 that are in array2
SweepParam = SweepParam(LogicalIndex);

for Index=length(SweepParam):-1:1

    load(FileList(Index).name, 'Field');

    valid=(Field.time>=TStart & Field.time<=TEnd).*(rem(1:length(Field.time), Subsample)==0)==1;
    t=Field.time(valid);
    y=Field.EL(valid);
    Y=fftshift(fft(y))/length(y);
    P(:, Index)=10*log10(abs(Y).^2);

    dt=t(2)-t(1);
    f(:, Index)=linspace(-0.5/dt,+0.5/dt,length(t));

end

%choose the name of the file
% SaveString=sprintf('OS_LEF1.2-DFB9-SHBYES-Rfront0.32-Rback0.99.mat');
% save(SaveString, 'f', 'P', '-mat');

%% Figure OS Map 

if isMap
    figPSM=figure();
    
    % Resize the matrix and extract central elements because 
    % otherwise the computational time is infinite:
    % I take 1 element every 64
    FreqMin = -1200; %GHz
    FreqMax = 1200; %GHz
    IndexMin=find(f(:,1)>=FreqMin, 1);
    IndexMax=find(f(:,1)>=FreqMax, 1);
    resizedf = f(IndexMin:64:IndexMax, :);
    resizedP = P(IndexMin:64:IndexMax, :);
    
    % Define the custom colormap
    % colorRange = [0  0 0.2430  % Blue
    %         0 1 1;   % Cyan
    %         0 1 0;   % Green
    %         1 1 0;   % Yellow
    %         1 0.6 0.1;   % Orange
    %         1 0 0];  % Red  % Color range corresponding to the output power range
    % 
    % % Define the power thresholds for color interpolation
    % thresholds = [max(max(P))-55; 
    %                 max(max(P))-50;
    %                 max(max(P))-45;
    %                 max(max(P))-30;
    %                 max(max(P))-15;
    %                 max(max(P))]; 
    % 
    % % Create a color axis using the custom colormap and range of values
    % customColormap=interp1(thresholds, colorRange, linspace(min(thresholds), max(thresholds), 256));
    
    [~,Ymesh] = meshgrid(resizedP(:,1)',SweepParam');
    contourf(resizedf',Ymesh,resizedP',200,'edgecolor','none');
    shading interp;
    % Create a color axis using the custom colormap and range of values
    %colormap(customColormap);
    %clim([-75 max(max(resizedP))])
    colormap('turbo');
    set(gca,'YDir','normal'); % to display the Y-axis in ascending order
    xlabel('Frequency (GHz)','fontsize',30,'interpreter','latex');
    
    if strcmp(SweepParameter, 'Current')==1
            ylabel('Current (mA)','fontsize',30,'interpreter','latex');
    elseif strcmp(SweepParameter, 'LEF')==1
            ylabel('LEF','fontsize',30,'interpreter','latex');
    elseif strcmp(SweepParameter, 'kDFB')==1
            ylabel('$\kappa_{DFB}$','fontsize',30,'interpreter','latex');
    end
    
    
    %xlim([-1000 1000])
    %ylim([280 600])
    
    set(gca,'xminortick', 'on', ...
        'yminortick', 'on','FontSize', 30,... 
            'TickLabelInterpreter','latex', 'TickDir', 'in') %Ticks
    
    c=colorbar;
    set(c,'TickLabelInterpreter', 'latex');
    c.Label.String='Power (dB)';
    c.Label.Interpreter='latex';
    c.Label.FontSize=30;
    
    set(gcf, 'color', 'w') ;
    set(figPSM,'GraphicsSmoothing','on')
    figPSM.Position=[(ScreenSize(1)-Width)/2 (ScreenSize(2)-Height)/2 Width Height];

end



end