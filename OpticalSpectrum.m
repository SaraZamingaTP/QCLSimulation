%% Representation of optical spectrum 

function OpticalSpectrum(TStart, TEnd, alphaf, kcoup)

%% Parameters setting

TStartfile=0;

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

%%
figOS=figure();

for current_index=length(InputCurrent):-1:1
    
    filenamestr=sprintf('Res_%gmAFrom0nsTo%gns_alphaf%g_kcoup%g.mat', ...
        InputCurrent(current_index), TEnd, alphaf, kcoup);
    filename = fullfile(pwd, filenamestr);
    load(filename);

    valid=(Field.time>=TStart & Field.time<=TEnd).*(rem(1:length(Field.time), Subsample)==0)==1;
    t=Field.time(valid);
    dt=t(2)-t(1);
     
    f=linspace(-0.5/dt,+0.5/dt,length(t));

    y=Field.EL(valid);
    Y=fftshift(fft(y))/length(y);
    
    P=abs(Y).^2;

    pl(current_index) = plot(f,P, 'linewidth', 1.5, ...
        'DisplayName', [num2str(InputCurrent(current_index)), 'mA']);
   
    hold on;

end

%axis settings
%xlim([-9.344 -9.340])
ylim([0 10])
xlabel('Frequency [GHz]','fontsize',30,'interpreter','latex')
ylabel('Optical Spectrum [dB]','fontsize',30,'interpreter','latex')
set(gca,'xminortick', 'on','yminortick', 'on','FontSize', 25,...
        'TickLabelInterpreter','latex', 'TickDir', 'in') %Ticks
set(gca,'linewidth', 1.5 , 'FontSmoothing', 'on', 'fontweight', 'normal')

%legend settings
legend(pl, 'location', 'northeast','NumColumns', 1,'fontsize',25,'interpreter','latex')
box off
legend boxoff


% Figure settings
ScreenSize=get(groot,'ScreenSize');
ScreenSize=ScreenSize(3:4);
Height=600;
Width=1200;
set(gcf, 'color', 'w') ;
set(figOS,'GraphicsSmoothing','on')
figOS.Position=[(ScreenSize(1)-Width)/2 (ScreenSize(2)-Height)/2 Width Height];

end
