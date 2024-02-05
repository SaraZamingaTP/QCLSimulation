% function to save in a .mat file the Optical Spectra: f, P, and Parameters
% used for Sweeping
function OpticalSpectra(TStart, TEnd)

Subsample=1;
FilenamePrefix = 'Res_'; 
% use the dir function to find all files that match the pattern
FileList = dir(fullfile(pwd, [FilenamePrefix, '*']));

%for loop to save the sweeping vector
for Index=length(FileList):-1:1
    BiasCurrent(Index,1)=str2double(extractBetween(FileList(Index).name,'_','mA'));
    LEF(Index,1)=str2double(extractBetween(FileList(Index).name,'alphaf','_'));
    kDFB(Index,1)=str2double(extractBetween(FileList(Index).name,'kcoup','_'));

    load(FileList(Index).name, 'Field');

    valid=(Field.time>=TStart & Field.time<=TEnd).*(rem(1:length(Field.time), Subsample)==0)==1;
    t=Field.time(valid);
    y=Field.EL(valid);
    Y=fftshift(fft(y))/length(y);
    P(:, Index)=10*log10(abs(Y).^2);

    dt=t(2)-t(1);
    f(:, Index)=linspace(-0.5/dt,+0.5/dt,length(t));

end

FreqMin = -1000; %GHz
FreqMax = 1000; %GHz
IndexMin=find(f(:,1)>=FreqMin, 1);
IndexMax=find(f(:,1)>=FreqMax, 1);
%resized f and P
Frequency = f(IndexMin:IndexMax, :);
PowerdB = P(IndexMin:IndexMax, :);

% Specify the directory where you want to save the figures
%SaveDir = 'D:\PhD\Simulation\GitSandbox\Results\Data\LICurve';
cd D:\PhD\Simulation\GitSandbox\Results\Data;
%choose the name of the file
SaveString=sprintf('OS_LEF1.2-DFB9-SHBYES-Rfront0.05-Rback0.99.mat');
save(SaveString, 'BiasCurrent', 'LEF', 'kDFB', 'Frequency', 'PowerdB', '-v7.3', '-mat');


end