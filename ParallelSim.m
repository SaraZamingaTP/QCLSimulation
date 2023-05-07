%launch TDTW_rungeFP.m for multiple simulations

function SimState=ParallelSim(isAlphaVar, isKcouplingVar)

%% input parameters for the simulation
%change these parameters for simulating different input conditions:
Folder=pwd; %where the files are 
TStart=0; % start time of the simulation [ns]
TEnd=200; %stop time [ns]
LoadStateString=''; %in case a new simulation should start from the final 
% state of a previous simulation, load here the relative Status file
% Define range of input values
inputCurr = 525:25:800; 

%% Simulation 

try

    if isAlphaVar==false && isKcouplingVar==false
        parfor i=1:length(inputCurr)
            ResTD{i}=Main_TDTW_DFBQCL(Folder, inputCurr(i),TStart,TEnd, LoadStateString, ...
                isAlphaVar, isKcouplingVar);
        end
    else

        for i=1:length(inputCurr)
            ResTD{i}=Main_TDTW_DFBQCL(Folder, inputCurr(i),TStart,TEnd, LoadStateString, ...
                isAlphaVar, isKcouplingVar);
        end
    end
    
%     FilenamePrefix = 'Res_';
% 
%     % use the dir function to find all files that match the pattern
%     FileList = dir(fullfile(pwd, [FilenamePrefix, '*']));
%     for i=numel(FileList):-1:1
%         filename = fullfile(pwd, FileList(i).name);
%         load(filename);
%         post_proc_simple_NEW(Field, Carriers, 2, 100);
%     end

    disp('Simulations completed\n!');
    SimState=true;
    
catch ME
    disp('Error occurred\n!');
    disp(ME.message);

    SimState=false;
end 
    

end
