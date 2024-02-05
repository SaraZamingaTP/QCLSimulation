
function AutocorrelationComputation(filename)
    %% Parameters of simulation
    Subsample=1;
    V=2240;
    L=2000; 
    Rfront=0.32;
    neff=3.3;
    TStart=1; 
    TEnd=500;

    load(filename)
    
    %% Generation of the plot of Power vs time
    valid=(Field.time>=TStart & Field.time<=TEnd).*(rem(1:length(Field.time), Subsample)==0)==1;
    t=Field.time(valid); %ns
    y=Field.EL(valid);
    PLtime=(0.5*2.997925e5*neff*8.854e-18/(1.6e-19)*...
                (abs(y).^2)*(1-Rfront)*(1/0.127)*V/L)*(1.6e-7);
    
    dt=t(2)-t(1); %around 30 fs, but unit of measure of the result is ns
     
    Limit=length(PLtime);
    
    MeanP = mean(PLtime);
    %shift with respect to the mean value
    P = PLtime - MeanP; %mW
    StdP = std(P);
    % Required to derived the auto-correlation of master 
    [intP,lagP] = xcorr(P,P);
    % lagP = lag index
    lagP=dt*lagP; %from index to timescale
    intP = intP/(StdP*StdP*Limit);
    
    %% plot autocorrelation
    figCorr=figure();
    plot(lagP,intP, ...
        'Linewidth', 2, 'color', [0.09,0.11,0.52]);
    
    xlabel('Time lag (ns)','fontsize',30,'interpreter','latex')
    ylabel('Norm. Intensity (arb. u.)','fontsize',30,'interpreter','latex')
    title('Auto Correlation','fontsize',30,'interpreter','latex')
    %xlim([-20 20])
    set(gca,'xminortick', 'on','yminortick', 'on', ...
    'FontSize', 30, 'TickLabelInterpreter','latex', 'TickDir', 'in'); %Ticks
    set(gcf, 'color', 'w') ;
    set(figCorr,'GraphicsSmoothing','on')
    
    %% Figure settings
    ScreenSize=get(groot,'ScreenSize');
    ScreenSize=ScreenSize(3:4);
    Height=600;
    Width=900;
    figCorr.Position=[(ScreenSize(1)-Width)/2 (ScreenSize(2)-Height)/2 Width Height];

end