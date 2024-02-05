function dati=DeviceData_OptExpr_FP() %#codegen
NumQDPopulations=1;            %[1,1] Number of sublevels populations (Note: this number must be odd!!!) []
%EnergyGapWL=1.1539;             %[1,1] Energy gap between CB and VB band edges in the WL [eV]
EnergyGapQDCentralPop=0.9879;   %[1,1] Characteristic interband transition energies for the QD central (most likely) QD population [eV]
InhomDeltaE=0.040/sqrt(8*log(2));
EnergyGapQD=[...                %[1,NumPops] Characteristic interband transition energies for the QD populations [eV]
    linspace(EnergyGapQDCentralPop+3*InhomDeltaE,EnergyGapQDCentralPop-3*InhomDeltaE,NumQDPopulations)];

EMass=9.11e-31;                 %[1,1] Electron mass [kg]
%% Device parameters
dati=struct(...
'L',            2000,...        %[1,1] total cavity length [um]
'nr',           3.3,...         %[1,1] effective refractive index []
'Func_alfa_i',  @Alfa_i,...     %[f(z)->1,1] Intrinsic waveguide losses [cm^-1]
'T',            273+20,...      %[1,1] temperature [K]
'GammaY',       @GammaY,...     %[f(z)->1,NumSlice] ridge transverse confinement factor []
'gamma_x',      0.127,...       %[1,1] field confinement factor in the QCL layers in the growth direction
'Func_Elettrodo',@Elettrodo,... %[f(z)->1,NumSlice integer] Identification of electrode names []
'NumQDPopulations',NumQDPopulations,...     %[1,1] Number of sublevels population (Note: this number must be odd!!!) []
'Gain',         2*[0.1363*0.77*1.05*0.8;0.1363/2*0.77*1.05]*1e-15,... % moltiplicazione per 0.7 in ES1 per tenere conto degli integrali di sovrapposizione
'Func_HomogeneousBroadeningGS',@Func_HomogeneousBroadeningGS,...               %[1,1] 
'tau_d',        1e-4,... %(ns)
'taue',         1e-3,... %(ns)
'R0',           0.95,...%.32,...         %[1,1] power reflectivity in z=0 facet
'RL',           0.3,...%.32,...         %[1,1] power reflectivity in z=L facet
'Np',           50,...
'omega_0',      2*pi*3*10^8/(9.34*10^(-6))/10^(9),... %(GHz) Reference frequency
'alphaf',       0.8,...
'f0',           1.1*10^(-7),...   %um^3
'Vol',          2240,...  %(um^3)
'Gamma',         0.15,... %(Ghz) Half Gain Linewidth %'Gamma',        320,... %(Ghz) Half Gain Linewidth
'epsb',         3.3^2);
end

function d=Width(z)
    d=ones(size(z))*5;
end

function A=Alfa_i(z)
    A=ones(size(z))*3.8;             
end

%Ntoth: numero di portatori totali in ogni slice
% function H=Func_HomogeneousBroadeningGS(Ntoth)
% H=7.60000000000000e+003.*ones(size(Ntoth));
% end

function Gy=GammaY(z)
    Gy=ones(size(z))*0.95;
end

function E=Elettrodo(z)
    E=ones(size(z));
end