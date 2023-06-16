function PowerSpectrumMap(TStart, TEnd, alphaf, kcoup, etai)

%Laser spectrogram as a function of the current 

%% % Load your data into Matlab
Subsample=1;
FilenamePrefix = 'Res_'; 
% use the dir function to find all files that match the pattern
FileList = dir(fullfile(pwd, [FilenamePrefix, '*']));

%for loop to save InputCurrent vector
for CurrentIndex=length(FileList):-1:1
    VectorParam=str2double(extractBetween(FileList(CurrentIndex).name,'_','mA'));
    InputCurrent(CurrentIndex)=VectorParam(1);
end

%sort in ascendent order the InputCurrent values
InputCurrent=sort(InputCurrent);

for CurrentIndex=length(InputCurrent):-1:1

    Filenamestr = sprintf('Res_%gmAFrom0nsTo%gns_alphaf%g_kcoup%.1f_etai0.8.mat', ...
        InputCurrent(CurrentIndex), TEnd, alphaf, kcoup); 
    filename = fullfile(pwd, Filenamestr);
    load(filename);

    valid=(Field.time>=TStart & Field.time<=TEnd).*(rem(1:length(Field.time), Subsample)==0)==1;
    t=Field.time(valid);
    y=Field.EL(valid);
    Y=fftshift(fft(y))/length(y);
    P(:, CurrentIndex)=10*log10(abs(Y).^2);

    %define just once the time and frequency vector, since it is in common
    %to all the current values
    if CurrentIndex==1
        dt=t(2)-t(1);
     % Define the sampling frequency
        %fs=1/dt;
        f(CurrentIndex,:)=linspace(-0.5/dt,+0.5/dt,length(t));
    end

end

%% Figure 
figPSM=figure();

% Define a custom colormap with 6 colors
cmap = [0.3 0.1 0.6  % Blue
        0 1 1;   % Cyan
        0 1 0;   % Green
        1 1 0;   % Yellow
        1 0.6 0.1;   % Orange
        1 0 0];  % Red
    
% Set the range of values for each color in the colormap
caxis_values = [max(max(P))-22; 
                max(max(P))-18;
                max(max(P))-17;
                max(max(P))-16;
                max(max(P))-15;
                max(max(P))]; 

% Create a color axis using the custom colormap and range of values
colormap(interp1(caxis_values, cmap, linspace(min(caxis_values), max(caxis_values), 256)));
% Set the colormap and scaling
clims=[max(max(P))-45, max(max(P))];
InputCurrent=InputCurrent';

% % Reshape the frequency and power spectrum matrices
% freq_data_2D = reshape(f(1,:), size(f, 1), []);
% for CurrentIndex=length(InputCurrent):-1:1
%     power_spectrum_db_2D = reshape(P(:, CurrentIndex), ...
%         size(P(:, CurrentIndex), 1), []);
% end
% 
% % Create a power spectrum map using surf
% surf(InputCurrent, freq_data_2D, power_spectrum_db_2D, clims);

imagesc(InputCurrent, f(1,:)', P, clims);
shading interp;
set(gca,'YDir','normal'); % to display the Y-axis in ascending order
xlabel('Bias Current (mA)','fontsize',30,'interpreter','latex');
ylabel('Frequency (GHz)','fontsize',30,'interpreter','latex');
ylim([-200 30])

set(gca,'xtick', InputCurrent,'xminortick', 'on', ...
    'yminortick', 'on','FontSize', 25,...
        'TickLabelInterpreter','latex', 'TickDir', 'in') %Ticks

c=colorbar;
set(c,'TickLabelInterpreter', 'latex');
c.Label.String='Power (dB)';
c.Label.Interpreter='latex';
c.Label.FontSize=30;


% Figure settings
ScreenSize=get(groot,'ScreenSize');
ScreenSize=ScreenSize(3:4);
Height=600;
Width=1200;

set(gcf, 'color', 'w') ;
set(figPSM,'GraphicsSmoothing','on')
figPSM.Position=[(ScreenSize(1)-Width)/2 (ScreenSize(2)-Height)/2 Width Height];


end