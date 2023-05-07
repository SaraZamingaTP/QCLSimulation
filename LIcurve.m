%% Function to represent the LIV curve 
function LIcurve(TStart, TEnd)

Subsample=1;
V=2400;
L=2000;

%% Computation of Output Power for the current sweep

FilenamePrefix = 'Res_'; 
% use the dir function to find all files that match the pattern
FileList = dir(fullfile(pwd, [FilenamePrefix, '*']));

for CurrentIndex=length(FileList):-1:1

    VectorParam=str2double(extractBetween(FileList(CurrentIndex).name,'_','mA'));
    InputCurrent(CurrentIndex)=VectorParam(1);
end

InputCurrent=sort(InputCurrent);
for CurrentIndex=length(InputCurrent):-1:1

    Filenamestr = sprintf('Res_%gmAFrom0nsTo300ns_alphaf0.8_kcoup6.5.mat', InputCurrent(CurrentIndex)); 
    filename = fullfile(pwd, Filenamestr);
    load(filename);
    valid=(Field.time>=TStart & Field.time<=TEnd).*(rem(1:length(Field.time), Subsample)==0)==1;
    t=Field.time(valid);
    y=Field.EL(valid);
    
    PLtime=(0.5*2.997925e5*3.3*8.854e-18/(1.6e-19)*...
            (abs(y).^2)*(1-0.3)*(1/0.127)*V/L)*(1.6e-8);
            
    %PLtime=0.5*0.127*2.997925e5*3.6*8.854e-18/(1.6e-19)*(1-0.3)*14*(2/25)*abs(y).^2*1.6*10^(-7);%.*1.6e-19*1/ResTD.Sim.dt/1e-9.*1e3;
    PLtimeAv(CurrentIndex)=mean(PLtime(t>t(end)-1));

end

%% Generation of the plot of Power vs Current
figLI=figure();
plot(InputCurrent,PLtimeAv, 'linewidth', 1.5, 'color', [0.191 0.18 0.56]); 

%axis settings
xlabel('Current [mA]','fontsize',30,'interpreter','latex');
ylabel('Output power [mW]','fontsize',30,'interpreter','latex');

set(gca,'xtick', InputCurrent,'xminortick', 'on','yminortick', 'on','FontSize', 25,...
        'TickLabelInterpreter','latex', 'TickDir', 'in') %Ticks
% Figure settings
ScreenSize=get(groot,'ScreenSize');
ScreenSize=ScreenSize(3:4);
Height=600;
Width=1200;

set(gcf, 'color', 'w') ;
set(figLI,'GraphicsSmoothing','on')
figLI.Position=[(ScreenSize(1)-Width)/2 (ScreenSize(2)-Height)/2 Width Height];


end