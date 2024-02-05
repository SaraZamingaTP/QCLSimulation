%launch TDTW_rungeFP.m for multiple simulations

function SimState=ParallelSim(isCurrentSweep, isAlphaVar, isKcouplingVar, isEtaiVar)

    isValidSimulation=true;

   %% check on the input parameters
    if nargin == 0  %in case no parameters are provided, the current sweep mode is activated
        isAlphaVar=false;
        isKcouplingVar=false;
        isEtaiVar=false;
        isCurrentSweep=true; 
        inputCurr=[310:1:314]; %default range for current bias sweep
        
    elseif nargin == 1 %if one param is provided, it is isCurrentSweep

        isAlphaVar=false;
        isKcouplingVar=false;
        isEtaiVar=false;

        if isCurrentSweep==false
            inputCurr = 300;
        else
            inputCurr= 250:25:800; %default range for current bias sweep
        end

    elseif nargin == 2

        isKcouplingVar=false;
        isEtaiVar=false;

        if (isCurrentSweep==true && isAlphaVar==true) || ...
            (isCurrentSweep==false && isAlphaVar==false)

            disp('Invalid input arguments\n');
            isValidSimulation=false;
        
        else
            
            if isCurrentSweep==true
                inputCurr= 250:25:800;
            else
                inputCurr=328; %default value
            end
            
        end

    elseif nargin == 3

        isEtaiVar=false;

        if (isCurrentSweep==true && (isAlphaVar==true || isKcouplingVar==true)) || ...
            (isCurrentSweep==false && isAlphaVar==false && isKcouplingVar==false) 

            disp('Invalid input arguments\n');
            isValidSimulation=false;

        else 

            if isCurrentSweep==true
               inputCurr= 250:25:800;
            else
                inputCurr=328; %default value
            end
        end

    elseif nargin == 4

        if (isCurrentSweep==true && (isAlphaVar==true ...
                || isKcouplingVar==true || isEtaiVar==true)) || ...
            (isCurrentSweep==false && isAlphaVar==false && ...
            isKcouplingVar==false && isEtaiVar==false) 

            disp('Invalid input arguments\n');
            isValidSimulation=false;
        
        else 
            
            if isCurrentSweep==true
               inputCurr= 250:25:800;
            elseif isEtaiVar==true
                inputCurr=400; %default value
            end
        end

    elseif nargin > 4

        disp('Too many input arguments\n');
        isValidSimulation=false;
    end

%% input parameters for the simulation
%change these parameters for simulating different input conditions:
Folder=pwd; %where the files are 
TStart=0; % start time of the simulation [ns]
TEnd=230; %stop time [ns]

LoadStateString=''; %in case a new simulation should start from the final 
% state of a previous simulation, load here the relative Status file
% Define range of input values


%% Simulation 

if isValidSimulation %only if the number of valid parameters is the good one
    
    try
    
        if isCurrentSweep %current sweep
            parfor i=1:length(inputCurr)
                ResTD{i}=Main_TDTW_DFBQCL(Folder, inputCurr(i),TStart,TEnd, LoadStateString, ...
                    isAlphaVar, isKcouplingVar, isEtaiVar);
            end

        else
    
            for i=1:length(inputCurr)
                ResTD{i}=Main_TDTW_DFBQCL(Folder, inputCurr(i),TStart,TEnd, LoadStateString, ...
                    isAlphaVar, isKcouplingVar, isEtaiVar);
            end
        end
    
        disp('Simulations completed\n!');
        SimState=true;
        
    catch ME
        disp('Error occurred\n!');
        disp(ME.message);
    
        SimState=false;
    end 
end

end
