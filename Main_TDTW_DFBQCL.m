function [ResTD, isSimComplete]=Main_TDTW_DFBQCL(Folder, Current, TStart, TEnd, ...
    LoadStateString, isAlphaVar, isKcouplingVar)
cd (Folder);
%%

% TStart, TEnd = starting and ending time of simulation in ns

if TEnd-TStart<=0
    throw  ('The simulation duration required is <=0ns!');
end

ResString1=sprintf('Res_%gmAFrom%gnsTo%gns',Current,TStart,TEnd);  
if exist(ResString1,'file')
    error ('File %s already exists!',ResString1);
end

SaveStateString1=sprintf('Status_%gmA_%gns',Current,TEnd);              % name of the file for saving data
if exist(SaveStateString1,'file')
    error ('File %s already exists!',SaveStateString1);
end

if isAlphaVar==false && isKcouplingVar==false

    alphaf=1.5; %alpha factor  
    kcoupling=5; %coupling coefficient

    SaveStateString2=sprintf("_alphaf%.1f_kcoup%.1f.mat", ...
                alphaf, kcoupling);
    SaveStateString = append(SaveStateString1,SaveStateString2);
    
    Dati=DeviceData_OptExpr_DFB(alphaf, kcoupling);                                                   % Dati is a structure with all the physical parameters
    
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
    
     ResTD=TDTW_runge_DFBQCL(Dati,Sim);
     ResString2=sprintf("_alphaf%.1f_kcoup%.1f.mat", alphaf, kcoupling);
     ResString = append(ResString1,ResString2);
     save(ResString, '-struct','ResTD');

     isSimComplete=true;

else
    if isAlphaVar==true

        alphaf=(0.8:0.2:1.2)'; %alpha factor vector 
        kcoupling=6.5; %coupling coefficient
    
        parfor alpha_index=1:length(alphaf)
    
            SaveStateString2=sprintf("_alphaf%.1f_kcoup%.1f.mat", ...
                alphaf(alpha_index), kcoupling);
            SaveStateString = append(SaveStateString1,SaveStateString2);
    
            Dati(alpha_index)=DeviceData_OptExpr_DFB(alphaf(alpha_index), kcoupling);                                                   % Dati is a structure with all the physical parameters
    
            Sim(alpha_index)= struct('dt',1*50e-6,'Length',TEnd-TStart,...                             % Sim is a structure with all the physical parameters
                'Pilot',@(e,t)Current*ones(size(e)),...
                'Pilot1',Current,...
                'Amplnoise_p',1.1*10^(-31),...
                'Amplnoise_r',1.1*10^(-31),...
                'PrintSimulationState',0.1,...
                'OutputPowerSampling',20,...
                'time_stop_Esp',10000,'OpticalPowerSpectrumSamplingStart',0,...
                'ComplexGainStep',1,'Mesh_factor',1,'reference_frequency','GS',...
                'subsampleRE',1,'LoadState',LoadStateString,'SaveState',SaveStateString);
            
            ResTD(alpha_index)=TDTW_runge_DFBQCL(Dati(alpha_index),Sim(alpha_index));
        
        end
    
        for alpha_index=1:length(ResTD)
            ResTD_NEW=ResTD(alpha_index);
            ResString2=sprintf("_alphaf%.1f_kcoup%.1f.mat", alphaf(alpha_index), kcoupling);
            ResString = append(ResString1,ResString2);
            save(ResString, '-struct','ResTD_NEW');

            isSimComplete=true;
        end
    
    end
    
    if isKcouplingVar==true
    
        alphaf=0.9; %alpha factor vector 
        kcoupling=(6:0.1:7)'; %coupling coefficient
       
        parfor k_index=1:length(kcoupling)
    
            SaveStateString2=sprintf("_alphaf%.1f_kcoup%.1f.mat", ...
                alphaf, kcoupling(k_index));
            SaveStateString = append(SaveStateString1,SaveStateString2);
            Dati(k_index)=DeviceData_OptExpr_DFB(alphaf, kcoupling(k_index));                                                   % Dati is a structure with all the physical parameters
    
            Sim(k_index)= struct('dt',1*50e-6,'Length',TEnd-TStart,...                             % Sim is a structure with all the physical parameters
                'Pilot',@(e,t)Current*ones(size(e)),...
                'Pilot1',Current,...
                'Amplnoise_p',1.1*10^(-31),...
                'Amplnoise_r',1.1*10^(-31),...
                'PrintSimulationState',0.1,...
                'OutputPowerSampling',20,...
                'time_stop_Esp',10000,'OpticalPowerSpectrumSamplingStart',0,...
                'ComplexGainStep',1,'Mesh_factor',1,'reference_frequency','GS',...
                'subsampleRE',1,'LoadState',LoadStateString,'SaveState',SaveStateString);
            
            ResTD(k_index)=TDTW_runge_DFBQCL(Dati(k_index),Sim(k_index));
        
        end
    
        for k_index=1:length(ResTD)
            ResTD_NEW=ResTD(k_index);
            ResString2=sprintf("_alphaf%.1f_kcoup%.1f.mat", alphaf, kcoupling(k_index));
            ResString = append(ResString1,ResString2);
            save(ResString, '-struct','ResTD_NEW');

            isSimComplete=true;
        end
    end
end

