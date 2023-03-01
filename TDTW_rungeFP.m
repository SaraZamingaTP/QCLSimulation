function [Res]=TDTW_rungeFP(Data,Sim) %#codegen
%Copyright 2011 Mattia Rossetti, Paolo Bardella, Mariangela Gioannini,
%Lorenzo Columbo, Carlo Silvestri
%Dynamic simulation of QCL using TDTW model
%Input:
%   Data:   the data structure with the definition of the device
%   Sim:    the structure with the simulation parameters
%   Measures: the structure with the information on the data to save
%             Optional; the default setting just saves output powers
%             !!!!NOT GIVEN HERE
%
%Output:
%   Res:    a structure with the saved data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%STEP 1 %%%%%%
%%%%%%%%%%%%%%  INITIALIZATION PROCEDURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check input parameters
if nargin<2
    error('requires 2 or 3 parameters');
end
% 

Sim=IntegrateSimParameters(Sim);                                            % if needed integrate sim structures field
if Sim.CheckParametersOnly==true
    Res=[];
    return
end

SimDuration=Sim.TEnd-Sim.TStart;

if Sim.SplitTime<SimDuration    
    TEnd=0;
    NParts=ceil(SimDuration/Sim.SplitTime);
    OriginalSaveState=Sim.SaveState;
    OriginalSplitTime=Sim.SplitTime;
    Sim.SplitTime=1e99;
    for k=1:NParts
        fprintf('This is part %d/%d of a %gns multi-part simulation\n',k, NParts,SimDuration);
        TStart=TEnd;
        if k==NParts
            TEnd=TStart+SimDuration-(NParts-1)*OriginalSplitTime;
        else
            TEnd=TStart+OriginalSplitTime;
        end
        Sim.TStart=TStart;
        Sim.TEnd=TEnd;        
        if TStart==0
            Sim.LoadState='';
        else
            Sim.LoadState=Res.Sim.SaveState;
        end
        Sim.SaveState=OriginalSaveState;
        Res=DDoSimulationOpt(Data, Sim, Measures);
    end
end

if Sim.CheckParametersOnly==true
    Res=[];
    return
end

%% Definition of additional commonly used phisical constants
Constants=struct('h',4.1357e-6,...                                          %Plank constant; [eV ns] (Electron volt per nanosecond)
                 'hcut',4.1357e-6/(2*pi),...;                               %Normalized Plank constant
                 'e',1.6e-19,...                                            %Electron charge [C]
                 'Kb',8.617343e-5,...                                       %Boltzmann constant [eV/'K]
                 'c',2.997925e5,...                                        %Vacuum speed of light [um/ns] 3e8m/s=3e8*1e6um/1e9ns=3e5um/ns
                 'ep0',8.854e-18/(1.6e-19));                                         %Vacuum permittivity [e/(um*V)]

             
disp(['Starting simulation on ' datestr(now)]);

                                       
%% GLOBAL REPLACEMENTS
%Some structures fields are here replaced with sorter names to increase
%overall readibility. Some of them are converted to int32 to allow for the
%automatic conversion by MATLAB coder


DeviceLength=Data.L;                                                        %[1,1] Total device length [um]
Temperature=Data.T;                                                         %[1,1] Device temperature [Kelvin]
KbTemp=Constants.Kb*Temperature;                                            %[1,1] Product of Boltzmann constant and temperature [eV]
NumQDPopulations=int32(Data.NumQDPopulations);                              %[1,1] Number of sublevels populations. Should be odd. [integer]
assert(rem(NumQDPopulations,2)==1);
onesNumQDPopulations=ones(NumQDPopulations,1);                              %[NumPop,1] Vector of 1's used for data expansion []
MeshFactor=int32(Sim.Mesh_factor);                                           %[1,1] Mesh factor from Javalojes accelerated model [integer]


GroupVelocity=Constants.c/Data.nr;                                          %[1,1] Group velocity [um/ns]
r0=sqrt(Data.R0);                                                           %[1,1] Field reflection coefficient at the facet in z=0 []
rL=sqrt(Data.RL);                                                           %[1,1] Field reflection coefficient at the facet in z=L []
t0=sqrt(1-Data.R0);                                                         %[1,1] Field transmission coefficient at the facet in z=0 []
tL=sqrt(1-Data.RL);                                                         %[1,1] Field transmission coefficient at the facet in z=L []



%% Calculation of the discretization space step and time step
%Fractional number of discretized slices, rounded to closest integer

Tau_dephasing=Data.tau_d  ;                                  %[1,1] dipoles dephasing time [ns]

dtnorm=Sim.dt/Tau_dephasing;

GroupVelocitynorm=GroupVelocity*Tau_dephasing/DeviceLength
  
DeviceLengthnorm=1
 
NumSlices=int32(round(1/(Sim.Mesh_factor*GroupVelocitynorm*dtnorm)));%[1,1]


%Recalculate the space step accordingly:
dz=double(DeviceLength)/double(NumSlices);                                  %[1,1] Spatial step [um]
Halfdz=dz/2;                                                                %[1,1] Half space step [um]
dznorm=double(DeviceLengthnorm)/double(NumSlices)           ;              %[1,1] Normalized space step

Datadt=dz/(GroupVelocity);                                                 %[1,1] Time step [ns]
Datadtnorm=dznorm/(GroupVelocitynorm);                                     %[1,1] Normalized  Time step 

%Print a note if the new simulation Timestep differs from the original one 
if Datadt~=Sim.dt
    fprintf('Simulation time step dt changed from %gns to %gns\n', Sim.dt, Datadt);
end
 
    fprintf('Simulation normalized time step is %g Tau_dephasing', Datadtnorm);

%Print also the total number of slices
fprintf('The %gum long device has been split in %d longitudinal slices\n', DeviceLength, NumSlices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sim.dt=Datadt;                                                              %[1,1] Simulation time step [ns]
NumStepsSteps=int32(floor(SimDuration/Sim.dt));                             %[1,1] Number of temporal steps
t=Sim.TStart+(1:1:double(NumStepsSteps))*Sim.dt;                            %[1,NumSteps] Integration times [ns]
tnorm=(Sim.TStart+(1:1:double(NumStepsSteps))*Sim.dt)/Tau_dephasing;                            %[1,NumSteps] Integration times [ns]


%ATTENZIONE!!!!*******************************
%Coordinates of the center of each slice along the longitudinal direction
z=dz/2:dz:DeviceLength-dz/2;                                                %[1,NumSlices] Position of slice centers [um]
zvet=0:dz:DeviceLength;
NumSlices=int32(length(z));                                                 %[1,1] Number of slices=length of the vector z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Alpha=Data.Func_alfa_i(z);                                                   %[1,NumSlices] waveguide losses [cm^-1]
%**************************************************************************
GammaX=Data.gamma_x;                                                        %[1,1] x-dir confinement factor []
MaterialLossesHalfed=Data.Func_alfa_i(z)/2/1e4;                             %[1,NumSlices]Field material losses, from cm^-1 [um^-1] 




%% Electric fields and random noise
% The field is expressed in physical units [eV/sec]. 
NumberOfSavedTemporalSteps=int32(max([2,MeshFactor+1]));                    %[1,1] Number of time steps to keep track of. [integer] %NOTE: WE USE HERE MESH FACTOR=1

%Initialization of the field

Sprog=AllocateComplex(NumberOfSavedTemporalSteps,NumSlices);          
Sregr=AllocateComplex(NumberOfSavedTemporalSteps,NumSlices); 
%Initialization of the carriers
RhocbGS0=AllocateDouble(NumQDPopulations,NumSlices);
RhocbGSp=AllocateComplex(NumQDPopulations,NumSlices);
RhocbGSm=AllocateComplex(NumQDPopulations,NumSlices);
%Initialization of the polarization envelopes
polp=AllocateComplex(NumQDPopulations,NumSlices);
polm=AllocateComplex(NumQDPopulations,NumSlices);


%% Field confinement factor
GammaY=Data.GammaY(z);
sqrtGammaY=sqrt(GammaY);
GammaYExpanded=repmat(GammaY,NumQDPopulations,1);
GammaXY=GammaX*GammaYExpanded;

GammaXY=GammaX*GammaY(1);      %!NOTE HERE GAMMAXY BECOMES A SCALAR
 

%% Initialization of indices before simulation                                                    
pp=NumberOfSavedTemporalSteps;                                              %[1,1] Index of the currently calculated time "t" []
pp2=NumberOfSavedTemporalSteps-1;                                           %[1,1] Index of the previously calculated time "t-dt"[]


%% Import of saved status
                                                            %[1,1] Initial value of the random generator []

if ~isempty(Sim.LoadState)   
        load(Sim.LoadState,'Sprog','Sregr','polp','polm',...
    'RhocbGS0','RhocbGSp');
end

 


%STEP 2 %%%%%%
%%%%%%%%%%%%%% Program main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NextTPrint=Sim.PrintSimulationState+Sim.TStart;                             %[1,1]Time at which the elapsed and remaining times are printed [ns]
tic;

NextTPrint1=0.01;
NextTPrint2=0.01;

Ep=AllocateComplex(NumQDPopulations,NumSlices);                                       %[NumQDPopulations,NumSlices] forward electric field inizialization
Em=AllocateComplex(NumQDPopulations,NumSlices);                                       %[NumQDPopulations,NumSlices] backward electric field inizialization

for it=1:NumStepsSteps

    
    if it==1 & strcmp(Sim.LoadState,'')==1
        
 
      
    Ep=Sim.Amplnoise_p*(rand(NumQDPopulations,NumSlices)+i*rand(NumQDPopulations,NumSlices));
    
    Em=Sim.Amplnoise_r*(rand(NumQDPopulations,NumSlices)+i*rand(NumQDPopulations,NumSlices));
    
    
    elseif  it==1 & strcmp(Sim.LoadState,'')==0
    
          
    Ep=repmat(Sprog(pp2,:),NumQDPopulations,1);
      
    Em=repmat(Sregr(pp2,:),NumQDPopulations,1); 
        
    
    end
    
    if it ~= 1 
        
          
    Ep=repmat(Sprog(pp2,:),NumQDPopulations,1);
     
    Em=repmat(Sregr(pp2,:),NumQDPopulations,1);
     
    end
    
    
    CurrentTime=t(it);                                                      %[1,1]Current step "real" time [ns]
    
  
    
    if CurrentTime>=NextTPrint1

        ElapsedHumanTime=toc;                                                       %[1,1]Elapsed time from the beginning of the simulation [s]
        RemainingHumanTime=ElapsedHumanTime/(CurrentTime-t(1))*(t(end)-CurrentTime);%[1,1]Remaining time before the end of the simulation[s]
        disp('time(ns)');
        disp([num2str(CurrentTime)]);
 
         disp('N0 (um^-3)');
         disp([num2str(RhocbGS0(1,end))]);

         disp('Power(mW)');
        disp([num2str(0.5*GammaXY*Constants.c*Data.nr*Constants.ep0*(1-Data.RL)*2*(14/25)*(abs(Sprog(pp2,end))^2+abs(Sregr(pp2,end))^2)*(1.6*10^(-7)))]);
                        
        NextTPrint1=NextTPrint1+0.1;%[1,1]Time at which the elapsed and remaining times are printed [ns]
    end
    
    if CurrentTime>=NextTPrint2
        
    save(Sim.SaveState,'Sprog','Sregr','polp','polm',...
    'RhocbGS0','RhocbGSp')
    
    NextTPrint2=NextTPrint2+1.0;
    
    end
    

   
%%%%%%% CARRIERS and POLARIZATIONS scaled RATE EQUATIONS%%%%%%%%%%%%%%%%%%%
        

        
        I=Sim.Pilot1;                                                          % pump current [mA]
       
       
        
        %% ADIMENSIONAL COEFFICIENTS IN THE EQUATIONS %%%%%%%%%%%%%%%
    
        CC=-1i*((Data.omega_0)*GammaXY*Data.Np)/(2*Constants.c*Constants.ep0*Data.nr);
        
        DD=(Data.Gamma*(1+1i*Data.alphaf)*(-1-1i*Data.alphaf)/Data.tau_d)*1i*Data.f0*Constants.ep0*Data.epsb;
        
        FF=Data.Gamma*(1+1i*Data.alphaf)/Data.tau_d;
        
%Second order Runge-Kutta algorithm is used to solve the ODE equations for each slice         
        
%Variation of the carriers variables       

drhocbGS0=10^(-12)*I/(Constants.e*Data.Vol)-RhocbGS0/Data.taue-1/(2*Constants.hcut)*imag(conj(Ep).*polp+conj(Em).*polm);
drhocbGSp=-RhocbGSp/Data.taue+1i/(4*Constants.hcut)*(conj(Em).*polp-Ep.*conj(polm));

%%%%% Variation of the microscopic polarizations %%%%%
drhopolp=-FF*polp-DD*(RhocbGS0.*Ep+RhocbGSp.*Em); 
drhopolm=-FF*polm-DD*(RhocbGS0.*Em+conj(RhocbGSp).*Ep); 

%Update of the values of carriers and polarization (first part of the Runge-Kutta algorithm)
        RhocbGS01=RhocbGS0+(Sim.subsampleRE/2)*Datadt*drhocbGS0;
        
        RhocbGSp1=RhocbGSp+(Sim.subsampleRE/2)*Datadt*drhocbGSp;
      
        polp1=polp+(Sim.subsampleRE/2)*Datadt*drhopolp;
    
        polm1=polm+(Sim.subsampleRE/2)*Datadt*drhopolm;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    
drhocbGS0=10^(-12)*I/(Constants.e*Data.Vol)-RhocbGS01./Data.taue-1/(2*Constants.hcut)*imag(conj(Ep).*polp1+conj(Em).*polm1);
drhocbGSp=-RhocbGSp1./Data.taue+1i/(4*Constants.hcut)*(conj(Em).*polp1-Ep.*conj(polm1));

        
      
    
        %%%%% Variation of the microscopic polarizations %%%%%
        
        drhopolp=-FF*polp1-DD*(RhocbGS01.*Ep+RhocbGSp1.*Em); 
        drhopolm=-FF*polm1-DD*(RhocbGS01.*Em+conj(RhocbGSp1).*Ep);
        
 %Update of the values of carriers and polarization (second part of the Runge-Kutta algorithm)
    
 
        RhocbGS0=RhocbGS0+Sim.subsampleRE*Datadt*drhocbGS0;

       RhocbGSp=RhocbGSp+Sim.subsampleRE*Datadt*drhocbGSp;
        
      
        polp=polp+Sim.subsampleRE*Datadt*drhopolp;
    
        polm=polm+Sim.subsampleRE*Datadt*drhopolm;
        
        
%%%%%%%% EQUAZIONI CAMPO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   %Field Equations 
    
    Sprog(pp,2:end)= Sprog(pp2,1:end-1)-dz*((MaterialLossesHalfed(1:end-1)).*...
    Sprog(pp2,1:end-1)-CC*polp(:,1:end-1));              
    
     Sregr(pp,1:end-1)=Sregr(pp2,2:end)-dz*((MaterialLossesHalfed(2:end)).*...
     Sregr(pp2,2:end)-CC.*polm(:,2:end));

    
    %%%%BOUNDARY CONDITIONS for a FP cavity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Sprog(pp,1)= r0*Sregr(pp,1)    ; 

    Sregr(pp,end)= rL*Sprog(pp,end);
    
 %   Boundary condition ring bidirectional
    
%    Sprog(pp,1)= r0*Sprog(pp,end);     

%    Sregr(pp,end)= rL*Sregr(pp,1);
      


%% saving results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
   
   
    %  If it is required to save the complex field in output from the device facets:
    if double(it)*Datadt>Sim.OpticalPowerSpectrumSamplingStart
        Res.Field.EL(it)=Sprog(pp,end).';
        
        Res.Field.PLphys(it)=0.5*Constants.c*Data.nr*Constants.ep0*(1-Data.RL)*14*(2/25)*abs((Sprog(pp,end).'))^2/(1.6*10^(-7))  ;   %power intracavity at the end facet [mW=10^-3*J/sec]
        
        
        Res.Field.E0(it)=Sprog(pp,1).';
        Res.Field.time(it)=CurrentTime;
        Res.Carriers.rhoGS0(:,it)=RhocbGS0(:,end);
        Res.Carriers.rhoGSp(:,it)=RhocbGSp(:,end);
        Res.Polariz.pp(:,it)=polp(:,end);
        Res.Polariz.pm(:,it)=polm(:,end);
       
    end
    
  
    Sprog(pp2,1:end)=Sprog(pp,1:end);
    Sregr(pp2,1:end)=Sregr(pp,1:end);

end %fine ciclo tempo
     
        
if ~isempty(Sim.SaveState)
    fprintf('Saving status to %s...', Sim.SaveState);          %%%% NOTE: CONSISTENT WITH LOAD STATE?
    %randstate=rand('state');
    %randstate=rng; %#ok    
    save(Sim.SaveState,'Sprog','Sregr','polm','polp',...
    'RhocbGS0','RhocbGSp');
        
    fprintf('done\n');    
end


disp(['Simulation terminated in ' num2str(cputime) 's. on ' datestr(now)]);

end

function y=AllocateComplex(r,c)
    y=complex(AllocateDouble(r,c));
end

function y=AllocateComplex3D(r,c,d)
    y=complex(AllocateDouble3D(r,c,d));
end

function y=AllocateDouble(r,c)
    y=AllocateVector(r*c);
    y=reshape(y,r,c);
end

function y=AllocateDouble3D(r,c,d)
    y=AllocateVector(r*c*d);
    y=reshape(y,r,c,d);
end

function y=AllocateVector(n)
    y=zeros(n,1);
end

function y=ReplicateDouble(val, r, c)
    y=AllocateDouble(r,c);
    for rc=1:r*c
        y(rc)=val;
    end
end



function Sim=IntegrateSimParameters(Sim)

DefaultSim={...
    'PrintSimulationState',0.1;...
    'OutputPowerSampling',10;...
    'time_stop_Esp',1e99;...
    'OpticalPowerSpectrumSamplingStart',0;...
    'ComplexGainStep',1;...
    'Mesh_factor',1;...    
    'subsampleRE',1;...    
    'LoadState','';...    
    'SaveState','Status';...    
    'ResFileName','Res';...    
    'TStart',0;...
    'SplitTime',1e99;...
    'CheckParametersOnly',false;...
    'LegacyFermiLevels',true;...
    'LegacySpontEmission',true};

for k=1:size(DefaultSim,1)
    if ~isfield(Sim,DefaultSim{k,1})
        Sim.(DefaultSim{k,1})=DefaultSim{k,2};
    end
end

if isfield(Sim,'Length')
    Sim.TEnd=Sim.TStart+Sim.Length;
    Sim=rmfield(Sim,'Length');
end
SimulationDuration=Sim.TEnd-Sim.TStart;

%Check for invalid situations 
if SimulationDuration<=0
    error ([mfilename ':InvalidDuration'], 'The simulation duration required is <=0ns!');
end

if Sim.SplitTime<SimulationDuration 
    if isempty(Sim.ResFileName)
        error ([mfilename ':CannotSplit'], 'To split the simulation, you must specify a destination filename in Sim.ResFileName');
    else
        [~,~, ext] = fileparts(Sim.ResFileName);
        if strcmpi(ext,'mat')==1
            error ([mfilename ':CannotSplit'], 'To split the simulation, you must specify a destination filename without .mat extension in Sim.ResFileName');
        end
    end
    
    if isempty(Sim.SaveState)
        error ([mfilename ':CannotSplit'], 'To split the simulation, you must specify a destination filename in Sim.SaveState');
    else
        [~,~, ext] = fileparts(Sim.SaveState);
        if strcmpi(ext,'mat')==1
            error ([mfilename ':CannotSplit'], 'To split the simulation, you must specify a destination filename without .mat extension in Sim.SaveState');
        end
    end        
else



if ~isempty(Sim.ResFileName) 
    [~,~, ext] = fileparts(Sim.ResFileName);
    if strcmpi(ext,'mat')==0
        Sim.ResFileName=sprintf('%s_%gmAFrom%gnsTo%gns.mat',Sim.ResFileName,Sim.Pilot(1,0), Sim.TStart, Sim.TEnd);
    end
end

if ~isempty(Sim.SaveState)
    [~,~, ext] = fileparts(Sim.SaveState);
    if strcmpi(ext,'mat')==0  
        Sim.SaveState=sprintf('%s_%gmAAt%gns.mat',Sim.SaveState,Sim.Pilot(1,0), Sim.TEnd);
    end
end
end

if exist(Sim.SaveState,'file')
    error ([mfilename ':InvalidSaveStateFile'], 'File %s already exists!',Sim.SaveState);
end

if exist(Sim.ResFileName,'file')
    error ([mfilename ':InvalidResFileName'], 'File %s already exists!',Sim.ResFileName);
end

if ~isempty(Sim.LoadState) && ~exist(Sim.LoadState,'file')
    error ([mfilename ':InvalidLoadStateFile'], 'File %s does not exists!',Sim.LoadState);
end

end

  
