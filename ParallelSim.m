%try to launch TDTW_rungeFP.m for different current values

%% input parameters for the simulation
%change these parameters for simulating different input conditions:
Folder=pwd; %where the files are 
TStart=0; % start time of the simulation [ns]
TEnd=300; %stop time [ns]
LoadStateString=''; %in case a new simulation should start from the final 
% state of a previous simulation, load here the relative Status file
% Define range of input values
inputCurr = 260:40:700; %range of input currents [mA] to run in parallel multiple simulations

remote_cluster=false;
%flag to select how to execute the parallel simulations:
%1) Batch: when using a remote cluster (server)
%2) parfor: for simulations on my computer 

%% Simulation 

if remote_cluster==true %SERVER

    %catch and handle errors to prevent the program from crashing
    %'try-catch statements are only useful when they surround the code that may generate errors
    try
        for i = length(inputCurr):-1:1
            job(i) = batch(@Main_TDTW, 1, ...
                {Folder, inputCurr(i), TStart, TEnd, LoadStateString});
        end
        %the decreasing index of the for loop eliminates the warning concerning the
        %changing size of job at each iteration, since Matlab knows in this way
        %which is gonna be the size of job -> no need for previous initialization
        
        % Retrieve the output from each job using for
        for i = 1:length(job)
          TEnd=300; %I need this because TEnd becomes Inf after first Sim
          %wait for the batch job to complete and retrieve the result
          wait(job(i), 'finished');
          delete(job(i));
          filename = sprintf('Res_%dmAFrom%dnsTo%dns.mat',inputCurr(i),TStart,TEnd);
          load(filename);
          post_proc_simple;
        end

        clear job
        disp('Simulations completed\n!');
    
    catch ME
    
        disp('Error occurred\n!');
        disp(ME.message);
    
    end

else %PERSONAL COMPUTER 
    try
        parfor i=1:length(inputCurr)
            ResTD{i}=Main_TDTW(Folder, inputCurr(i),TStart,TEnd, LoadStateString);
        end

        for i=length(inputCurr):-1:1
            TEnd=300; 
            filename = sprintf('Res_%dmAFrom%dnsTo%dns.mat',inputCurr(i),TStart,TEnd);
            load(filename);
            post_proc_simple;
        end
        
    catch ME
        disp('Error occurred\n!');
        disp(ME.message);
    end 
    
end 

