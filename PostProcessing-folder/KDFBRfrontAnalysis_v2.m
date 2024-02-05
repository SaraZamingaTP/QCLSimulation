
%Data
kDFB= [6, 9, 12, 15, 20, 25, 30, 35, 40];
IthAR=[380, 351, 329, 307, 290, 280, 275, 270, 269];
Ithcleaved=[306,295,286,278,271,266,263,262,261];

DeltaIsmAR=[160, 164, 8, 16, 10, 5, 1, 1, 0];
DeltaIsmcleaved=[64, 81, 94, 105, 113, 107, 78, 58, 45];


%% Fit results
kExtended=0:0.5:40;
RExtended=0.00:0.01:0.32;

DeltaIsmcleavedINTERP = interp1(kDFB, DeltaIsmcleaved, kExtended, 'pchip');
DeltaIsmARINTERP = interp1(kDFB, DeltaIsmAR, kExtended, 'pchip');

fig1 = figure(1);
%% 1D- CrossCorrelation and AutoCorrelation computation
subplot(1,2,1);
LegendString=sprintf('R$_{front}$=0.32, R$_{back}$=0.99\n');

p(1)=plot(kDFB, DeltaIsmcleaved, 'o', ...
    'DisplayName', LegendString, ...
     "MarkerEdgeColor", "#2722A0", 'MarkerFaceColor',"#2722A0");
hold on
plot(kExtended, DeltaIsmcleavedINTERP, ...
    'linewidth', 1.5, 'Color', "#2722A0");
hold on

LegendString=sprintf('R$_{front}$=0.05, R$_{back}$=0.99\n');
p(2)=plot(kDFB, DeltaIsmAR, 'o', 'DisplayName', LegendString, ...
    "MarkerEdgeColor", "#FA004C", 'MarkerFaceColor',"#FA004C");
hold on
plot(kExtended, DeltaIsmARINTERP, ...
    'linewidth', 1.5, 'Color', "#FA004C");

xlabel('$\kappa_{DFB}$ (cm$^{-1}$)','fontsize',20,'interpreter','latex')
ylabel('$\Delta$ I$_{sm}$ (mA)', 'fontsize',20,'interpreter','latex')
xlim([0 40])
hold on

set(gca,'xminortick', 'on', ...
    'yminortick', 'on','FontSize', 20,... 
        'TickLabelInterpreter','latex', 'TickDir', 'in') %Ticks

legend(p, 'location', 'northeast',...
    'NumColumns', 1,'fontsize',20,'interpreter','latex')

subplot(1,2,2);

LegendString=sprintf('R$_{front}$=0.32, R$_{back}$=0.99\n');
plot(kDFB, Ithcleaved, ...
    'linewidth', 1.5, 'DisplayName', LegendString,'Color', "#2722A0", ...
    "Marker","diamond", "MarkerEdgeColor", "#2722A0", 'MarkerFaceColor',"#2722A0");
hold on
LegendString=sprintf('R$_{front}$=0.05, R$_{back}$=0.99\n');
plot(kDFB, IthAR, ...
    'linewidth', 1.5, 'DisplayName', LegendString,'Color', "#FA004C", ...
    "Marker","diamond", "MarkerEdgeColor", "#FA004C", 'MarkerFaceColor',"#FA004C");
xlim([6 40])
xlabel('$\kappa_{DFB}$ (cm$^{-1}$)','fontsize',20,'interpreter','latex')
ylabel('I$_{th}$ (mA)', 'fontsize',20,'interpreter','latex')
hold on
legend('location', 'northeast',...
    'NumColumns', 1,'fontsize',20,'interpreter','latex')
set(gca,'xminortick', 'on', ...
    'yminortick', 'on','FontSize', 20,... 
        'TickLabelInterpreter','latex', 'TickDir', 'in') %Ticks

%Display the value of AC  
%fprintf('AC=%.2f\n', max(abs(int2)));
set(gcf, 'color', 'w') ;

% Size
ScreenSize=get(groot,'ScreenSize');
ScreenSize=ScreenSize(3:4);
Height=500;
Width=1200;
set(fig1,'GraphicsSmoothing','on')
fig1.Position=[(ScreenSize(1)-Width)/2 (ScreenSize(2)-Height)/2 Width Height];

%%

DeltaIsm=[160, 164, 8, 16, 10, 5, 1, 1, 0; ... % cleaved
    64, 81, 94, 105, 113, 107, 78, 58, 45]; %AR

ft = 'pchipinterp';
for RowIndex=size(DeltaIsm, 1):-1:1
    [fitresult, ~] = fit(kDFB', DeltaIsm(RowIndex,:)', ft, 'Normalize', 'on' );
    DeltaIsminterp(RowIndex,:)=fitresult(kExtended);
end

RfrontVector=[0.05, 0.32];
for ColIndex=size(DeltaIsminterp,2):-1:1
    [fitresult, ~] = fit(RfrontVector', DeltaIsminterp(:,ColIndex), ...
        ft, 'Normalize', 'on' );
    DeltaIsminterp2(:,ColIndex)=fitresult(RExtended);
end

fig3DkDFBvsRfront=figure();
[X,Y] = meshgrid(kExtended,RExtended);


subplot(1,2,1);
surf(X,Y, DeltaIsminterp2, 'EdgeColor', 'none', 'FaceAlpha', 1);
colormap turbo;
%grid minor
% hold on
% plot3(repmat(9, [28,1]), RExtended,  DeltaIsminterp2(:,18), ...
%     'Linewidth', 2, 'Color', [1, 1, 0]);
% hold on
% plot3(kExtended, repmat(0.32, [81,1]),  DeltaIsminterp2(end,:), ...
%     'Linewidth', 2, 'Color', [1, 1, 0]);

shading interp;
%shadowplot x;
%shadowplot y;
%Axis settings
zlabel('$\Delta I_{sm}$ (mA)','fontsize',24,'interpreter','latex');
xlabel('$\kappa_{DFB}$ ($cm^{-1}$)','fontsize',24,'interpreter','latex');
ylabel('$R_{front}$','fontsize',24,'interpreter','latex');
set(gca,'xminortick', 'on','yminortick', 'on','zminortick', 'on', ...
    'FontSize', 24, 'TickLabelInterpreter','latex', 'TickDir', 'in');%, ...
    %'XTick', kDFB, 'YTick', RfrontVector, 'ZTick', 0:30:160); %Ticks
% xlim([0, 42])
% ylim([0, 0.4])
% zlim([0, 150])
%colorbar settings
% c=colorbar;
% set(c,'TickLabelInterpreter', 'latex', 'Limits', [0, 164]);
% c.Label.String='mA';
% c.Label.Interpreter='latex';
% c.Label.FontSize=24;

%Figure settings
% ScreenSize=get(groot,'ScreenSize');
% ScreenSize=ScreenSize(3:4);
% Height=600;
% Width=1300;

% set(gcf, 'color', 'w') ;
% set(fig3DkDFBvsRfront,'GraphicsSmoothing','on')
% fig3DkDFBvsRfront.Position=[(ScreenSize(1)-Width)/2 (ScreenSize(2)-Height)/2 Width Height];


%% 

Ith=[306,295,286,278,271,266,263,262,261; ... % cleaved
    380, 351, 329, 307, 290, 280, 275, 270, 269]; %AR

ft = 'pchipinterp';
for RowIndex=size(Ith, 1):-1:1
    [fitresult, ~] = fit(kDFB', Ith(RowIndex,:)', ft, 'Normalize', 'on' );
    Ithinterp(RowIndex,:)=fitresult(kExtended);
end

RfrontVector=[0.05, 0.32];
for ColIndex=size(Ithinterp,2):-1:1
    [fitresult, ~] = fit(RfrontVector', Ithinterp(:,ColIndex), ...
        ft, 'Normalize', 'on' );
    Ithinterp2(:,ColIndex)=fitresult(RExtended);
end


[X,Y] = meshgrid(kExtended,RExtended);
%fig3DkDFBvsRfront=figure();
subplot(1,2,2);
surf(X,Y, Ithinterp2, 'EdgeColor', 'none', 'FaceAlpha', 1);
colormap turbo;
%grid minor
% hold on
% plot3(repmat(9, [28,1]), RExtended,  DeltaIsminterp2(:,18), ...
%     'Linewidth', 2, 'Color', [1, 1, 0]);
% hold on
% plot3(kExtended, repmat(0.32, [81,1]),  DeltaIsminterp2(end,:), ...
%     'Linewidth', 2, 'Color', [1, 1, 0]);

shading interp;
%shadowplot x;
%shadowplot y;
%Axis settings
zlabel('$I_{th}$ (mA)','fontsize',24,'interpreter','latex');
xlabel('$\kappa_{DFB}$ ($cm^{-1}$)','fontsize',24,'interpreter','latex');
ylabel('$R_{front}$','fontsize',24,'interpreter','latex');
set(gca,'xminortick', 'on','yminortick', 'on','zminortick', 'on', ...
    'FontSize', 24, 'TickLabelInterpreter','latex', 'TickDir', 'in');%, ...
    %'XTick', kDFB, 'YTick', RfrontVector, 'ZTick', 0:30:160); %Ticks
% xlim([0, 42])
% ylim([0, 0.4])
% zlim([0, 150])
%colorbar settings
% c=colorbar;
% set(c,'TickLabelInterpreter', 'latex', 'Limits', [260, 500]);
% c.Label.String='mA';
% c.Label.Interpreter='latex';
% c.Label.FontSize=24;

%Figure settings
ScreenSize=get(groot,'ScreenSize');
ScreenSize=ScreenSize(3:4);
Height=600;
Width=1300;

set(gcf, 'color', 'w') ;
set(fig3DkDFBvsRfront,'GraphicsSmoothing','on')
fig3DkDFBvsRfront.Position=[(ScreenSize(1)-Width)/2 (ScreenSize(2)-Height)/2 Width Height];
