function ResTD=Main_TDTW(Folder, Current, TStart, TEnd, LoadStateString)

cd (Folder);
%%

% TStart, TEnd = starting and ending time of simulation in ns

if TEnd-TStart<=0
    throw  ('The simulation duration required is <=0ns!');
end

ResString=sprintf('Res_%gmAFrom%gnsTo%gns.mat',Current,TStart,TEnd);  
if exist(ResString,'file')
    error ('File %s already exists!',ResString);
end

SaveStateString=sprintf('Status_%gmA_%gns.mat',Current,TEnd);              % name of the file for saving data
if exist(SaveStateString,'file')
    error ('File %s already exists!',SaveStateString);
end


Dati=DeviceData_OptExpr_FP;                                                   % Dati is a structure with all the physical parameters


Sim= struct('dt',1*50e-6,'Length',TEnd-TStart,...                             % Sim is a structure with all the physical parameters
    'Pilot',@(e,t)Current*ones(size(e)),...
    'Pilot1',Current,...
    'Amplnoise_p',1.1*10^(-31),...
    'Amplnoise_r',1.1*10^(-31),...
    'PrintSimulationState',0.1,...
    'OutputPowerSampling',20,...
    'time_stop_Esp',10000,'OpticalPowerSpectrumSamplingStart',0,...
    'ComplexGainStep',1,'Mesh_factor',1,'reference_frequency','GS',...
    'subsampleRE',1,'LoadState',LoadStateString,'SaveState',SaveStateString);

ResTD=TDTW_rungeFP(Dati,Sim);
save(ResString,'-struct','ResTD');



    
end

