%% Function to compute the LI curve and save the data in a .mat file 
function LIcurve(TStart, TEnd, Rfront)
%% Simulations Data
V=2240;
L=2000;

%% Computation of Output Power for the current sweep

Subsample=1;
FilenamePrefix = 'Res_'; 
% use the dir function to find all files that match the pattern
FileList = dir(fullfile(pwd, [FilenamePrefix, '*']));

for CurrentIndex=length(FileList):-1:1

    VectorParam=str2double(extractBetween(FileList(CurrentIndex).name,'_','mA'));
    InputCurrent(CurrentIndex)=VectorParam(1);
end

InputCurrent=sort(InputCurrent);

for CurrentIndex=length(InputCurrent):-1:1

    load(FileList(CurrentIndex).name, 'Field');

    valid=(Field.time>=TStart & Field.time<=TEnd).*(rem(1:length(Field.time), Subsample)==0)==1;
    t=Field.time(valid);
    y=Field.EL(valid);
    
    PLtime=(0.5*2.997925e5*3.3*8.854e-18/(1.6e-19)*...
            (abs(y).^2)*(1-Rfront)*(1/0.127)*V/L)*(1.6e-7);
            
    PLtimeAv(CurrentIndex)=mean(PLtime(t>t(end)-1));

end

% Specify the directory where you want to save the figures
cd D:\PhD\Simulation\GitSandbox\Results\Data;
%choose the name of the file
SaveString=sprintf('LI_LEF1.2-DFB9-SHBYES-Rfront0.05-Rback0.99.mat');
save(SaveString, 'InputCurrent', 'PLtimeAv', '-mat');

end