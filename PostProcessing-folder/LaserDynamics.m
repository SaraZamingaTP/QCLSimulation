function LaserDynamics(Filename)

load(Filename);
TurnOnPower=-20; %dB
MaxSMSR = 30; % Reference Side-Mode-Suppression-Ratio (SMSR)

NbrSpectra=length(BiasCurrent); 
 
%vector string defining the dynamic behavior of the laser, for the triplet
% (BiasCurrent, LEF, kDFB)
State=cell(NbrSpectra, 1);
ColorState=cell(NbrSpectra, 1);
for SpectraIndex=NbrSpectra:-1:1

    %% Find the laser mode separation

    PeakHeight = -90;
    % Find peaks higher than the threshold
    [peaks, peakLocations] = findpeaks(PowerdB(:, SpectraIndex), Frequency(:, SpectraIndex), ...
        'MinPeakHeight', PeakHeight, 'MinPeakDistance', 20);
    LaserModeSeparation=[];
    for FreqIndex = length(peakLocations):-1:2
        diff = abs(peakLocations(FreqIndex) - peakLocations(FreqIndex-1));
        LaserModeSeparation = [diff, LaserModeSeparation];
    end
    LaserModeSeparationAv= mean(LaserModeSeparation);


    %% Find the if the laser is OFF, Single-mode or Multimode

    % Find the maximum power and its index.
    [maxPower, ~] = max(PowerdB(:, SpectraIndex));
    
    [SideModesIntensity, SideModesFreq] = findpeaks(PowerdB(:, SpectraIndex), Frequency(:, SpectraIndex), ...
        'MinPeakHeight', maxPower - MaxSMSR, 'MinPeakDistance', 20);
       
    for ModesIndex=length(SideModesIntensity):-1:1
        if SideModesIntensity(ModesIndex) < -60
            SideModesFreq(ModesIndex)=0;
        end
    end

    % Use logical indexing to eliminate elements equal to 0
    SideModesFreq = SideModesFreq(SideModesFreq ~= 0);
    
    % Check if the laser is off or on.
    if (maxPower < TurnOnPower && length(SideModesFreq)==1) ...
            || isempty(SideModesFreq)
        State{SpectraIndex}='OFF'; 
        ColorState{SpectraIndex}='#120C5F'; %navy blue
    else
    
        
        if length(SideModesFreq) == 1
            State{SpectraIndex}='SM'; 
            ColorState{SpectraIndex}='#FF7800'; %orange
        elseif length(SideModesFreq) > 1
            
            SideModesSeparation = []; 
            % Iterate through the vector in reverse order
            for FreqIndex = length(SideModesFreq):-1:2
                diff = abs(SideModesFreq(FreqIndex) - SideModesFreq(FreqIndex-1));
                SideModesSeparation = [diff, SideModesSeparation];
            end
            
            SideModesSeparationAv=mean(SideModesSeparation);
            
            % if the distance is around 22 GHz
            if round(SideModesSeparationAv/22) == 1
                State{SpectraIndex}='D'; 
                ColorState{SpectraIndex}='#1EF525'; %green
            else
             % if the distance is a multiple of 22 GHz
                State{SpectraIndex}='H'; 
                ColorState{SpectraIndex}='#BA1EF5'; %purple
            end

            
        end
        
    end

end 

cd D:\PhD\Simulation\GitSandbox\Results\Data;
%choose the name of the file
SaveString=sprintf('OS_LEF1.2-DFB9-SHBYES-Rfront0.05-Rback0.99-StateVersion.mat');
save(SaveString, 'BiasCurrent', 'LEF', 'kDFB', 'Frequency', 'PowerdB', 'State', 'ColorState', '-mat');

end