function MapDynamics(Filename)

FileList = dir(fullfile(pwd, '*-StateVersion.mat'));
%load(Filename, 'BiasCurrent', 'LEF', 'kDFB', 'State');

%unique function takes just once the elements that repeat themselves in a
%vector
CurrentVal=unique(BiasCurrent); 
LEFVal=unique(LEF);
kDFBVal=unique(kDFB);

%% Create State 3D-Matrix for each (BiasCurrent, LEF, kDFB) combination

xSize=length(CurrentVal);
ySize=length(LEFVal);
zSize=length(kDFBVal);
%Initialize 3D Matrix
StateMatrix=cell(xSize, ySize, zSize);
for CurrentIndex=xSize:-1:1
    for LEFIndex=ySize:-1:1
        for kDFBIndex=zSize:-1:1
            MatchIndex = find(BiasCurrent == CurrentVal(CurrentIndex) & ...
                LEF == LEFVal(LEFIndex) & kDFB==kDFBVal(kDFBIndex));
            StateMatrix{CurrentIndex, LEFIndex, kDFBIndex}=State{MatchIndex}; 
        end
    end
end


%% Create a colormap for the states
uniqueStates = unique(StateMatrix, 'stable'); 
%'stable' is needed otherwise the strings are saved in alphabetic order
%(D, H, OFF, SM) and not is order of occurrence (OFF, SM, H, D)

% Define a custom colormap for your states (in this example, we use a simple colormap)
%customColors = [0 0 1; 0 1 0; 1 0 0; 1 1 0]; % Blue, Green, Red, Yellow
customColors ={'#120C5F'; '#FF7800'; '#BA1EF5'; '#1EF525'};
% Map unique string states to color indices
colorIndices = zeros(size(StateMatrix));
for i = 1:length(uniqueStates)
    colorIndices(ismember(StateMatrix, uniqueStates{i})) = i;
end

for ColorIndex=length(customColors):-1:1
    hex = customColors{ColorIndex}(2:end); % Remove '#' symbol
    r = hex2dec(hex(1:2)) / 255;
    g = hex2dec(hex(3:4)) / 255;
    b = hex2dec(hex(5:6)) / 255;
    rgb(ColorIndex, :) = [r; g; b];
end

for ColorIndex=length(colorIndices):-1:1
    RGBVector(ColorIndex,:)=rgb(colorIndices(ColorIndex),:);
end

% Scatter points in 3D space with colors based on states
figMap=figure();
scatter3(CurrentVal(:), LEFVal(:), kDFBVal(:), 50, RGBVector, 'filled');%, 50, rgb, 'filled');
%scatter3(CurrentVal(:), LEFVal(:), kDFBVal(:), 50, colorIndices(:), 'filled');

% Create a colormap
colormap(rgb);
%colorbar;
shading interp;
xlabel('Current (mA)','fontsize',30,'interpreter','latex');
ylabel('LEF','fontsize',30,'interpreter','latex');
zlabel('$\kappa_{DFB} (cm^{-1})$','fontsize',30,'interpreter','latex');

% Create a legend
%legend(stateUnique, 'Location', 'EastOutside');
set(gca,'xminortick', 'on', ...
        'yminortick', 'on', 'zminortick', 'on', 'FontSize', 30,... 
            'TickLabelInterpreter','latex', 'TickDir', 'in') %Ticks
    

% set(c,'TickLabelInterpreter', 'latex');
% c.Label.String='Power (dB)';
% c.Label.Interpreter='latex';
% c.Label.FontSize=30;

set(gcf, 'color', 'w') ;
set(figMap,'GraphicsSmoothing','on')

%% Figure settings
ScreenSize=get(groot,'ScreenSize');
ScreenSize=ScreenSize(3:4);
Height=600;
Width=1200;
figMap.Position=[(ScreenSize(1)-Width)/2 (ScreenSize(2)-Height)/2 Width Height];

end