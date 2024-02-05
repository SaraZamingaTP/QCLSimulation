%struct with Fields:
% Rfront, LEF, kDFB, Ith, DeltaIsm, Pout
clear all

SimStruct=struct('Rfront', [], 'LEF', [], 'kDFB', [], ...
    'Ith', [], 'DeltaI', [], 'Pout', []);
fieldNames = fieldnames(SimStruct);
% Count the number of fields
NbrField = numel(fieldNames);
VectorCompleteField=zeros(1,6);
StructRow=1; 
NbrElements=0;
ErrorFlag=0; 
%% Reading the txt file with all the infos
%load Simulation Results
cd 'D:\PhD\Simulation\GitSandbox\Results\June2023';
fileID=fopen('SummaryResults.txt');

while ~feof(fileID)
    
    %if VectorCompleteField == 6, I can start to read a new data set
    if sum(VectorCompleteField)==6 
            %before going on, check that the Fields with vector as
            %variables, have vectors with the same size

            %reinitialize VectorFlag and increase StructRow 
            VectorCompleteField=zeros(1,6);
            NbrElements=0;
            if ErrorFlag==0 %if the data relative to a precise Rfront
                %and LEF are not correct, I just subscribe them
                StructRow=StructRow+1;
            end
            ErrorFlag=0;
           
    end

    line = fgetl(fileID);
    % Check if the line contains a variable assignment
    if contains(line, '=')
        % Extract the variable name and values
        [varName, values] = strtok(line, '=');
        varName = strtrim(varName);
        
        MatchIndex = find(strcmp(varName, fieldNames)==1);
        if VectorCompleteField(MatchIndex) == 0
           VectorCompleteField(MatchIndex)=1;
        else %if it is already at 1, it means there is an error in the file
            disp('Error in the file! Field already complete\n');
        end

        if values(2)=='[' %for vectors
            values = strtrim(values(3:end-1));
        else                %for scalars
            values = strtrim(values(2:end));
        end
        
        % Convert the values to a numeric array
        values = str2double(strsplit(values, ';'));
        if strcmp(varName, 'LEF') || strcmp(varName, 'Rfront') %if should be a scalar 
            if numel(values)>1
                sprintf('%s must be a scalar\n', varName);
            end
        else %if it is a Field with vector
            if NbrElements == 0 %initialization
                NbrElements = numel(values);
            else
                if numel(values)~=NbrElements
                    ErrorFlag=1;
                    sprintf('%s has incoherent size: expected %d elements!\n', varName, NbrElements);
                end
            end
        end
        % Save the values into the structure
        SimStruct(StructRow).(varName) = values;
    end
end

% Close the file
fclose(fileID);

%% Fit results
ft = 'pchipinterp';
kExtended=7:0.5:40;
RExtended=0.05:0.01:0.32;
%% Represent data in 3D

kDFB=SimStruct.kDFB; %need this otherwise in meshgrid command I get an error

%construct RfrontVector
for RowIndex=length(SimStruct):-1:1
    RfrontVector(RowIndex)=SimStruct(RowIndex).Rfront;
    %Z is a matrix with DeltaI values
    Z(RowIndex,:)=SimStruct(RowIndex).DeltaI;
    [fitresult, ~] = fit(kDFB', Z(RowIndex,:)', ft, 'Normalize', 'on' );
    Zinterp(RowIndex,:)=fitresult(kExtended);
end

for ColIndex=size(Zinterp,2):-1:1
    [fitresult, ~] = fit(RfrontVector', Zinterp(:,ColIndex), ...
        ft, 'Normalize', 'on' );
    Zinterp2(:,ColIndex)=fitresult(RExtended);
end

[X,Y] = meshgrid(kExtended,RExtended);
fig3DkDFBvsRfront=figure();
surf(X,Y,Zinterp2, 'EdgeColor', 'none', 'FaceAlpha', 1);
colormap turbo;
%grid minor
hold on
plot3(repmat(12, [28,1]), RExtended, Zinterp2(:,11), ...
    'Linewidth', 2, 'Color', [1, 1, 0]);
hold on
plot3(kExtended, repmat(0.05, [67,1]), Zinterp2(1,:), ...
    'Linewidth', 2, 'Color', [1, 1, 0]);

shading interp;
%shadowplot x;
%shadowplot y;
%Axis settings
zlabel('$\Delta I_{sm}$ (mA)','fontsize',30,'interpreter','latex');
xlabel('$\kappa_{DFB}$ ($cm^{-1}$)','fontsize',30,'interpreter','latex');
ylabel('$R_{front}$','fontsize',30,'interpreter','latex');
set(gca,'xminortick', 'on','yminortick', 'on','zminortick', 'on', ...
    'FontSize', 30, 'TickLabelInterpreter','latex', 'TickDir', 'in', ...
    'XTick', kDFB, 'YTick', RfrontVector, 'ZTick', 0:30:160); %Ticks
xlim([0, 42])
ylim([0, 0.4])
zlim([0, 150])
%colorbar settings
c=colorbar;
set(c,'TickLabelInterpreter', 'latex', 'Limits', [0, 150]);
c.Label.String='mA';
c.Label.Interpreter='latex';
c.Label.FontSize=30;

%Figure settings
ScreenSize=get(groot,'ScreenSize');
ScreenSize=ScreenSize(3:4);
Height=600;
Width=1300;

set(gcf, 'color', 'w') ;
set(fig3DkDFBvsRfront,'GraphicsSmoothing','on')
fig3DkDFBvsRfront.Position=[(ScreenSize(1)-Width)/2 (ScreenSize(2)-Height)/2 Width Height];
