% efficiency
isEtaIVar=true;
isAVar=false;

alpha_i= 3.8; %cm^-1, waveguide loss
Lcav=0.2; %cm, cavity length

if isAVar && ~isEtaIVar
    A=linspace(0, 4, 32)'; %experimental value (Coldren, page 148)
    eta_i=0.8;
    GammaGth=A./Lcav+alpha_i;
    eta_d=eta_i.*(A./(GammaGth.*Lcav));
elseif isEtaIVar && ~isAVar
    A=1.2;
    eta_i=linspace(0.4, 0.9, 20)';
    GammaGth=A./Lcav+alpha_i;
    eta_d=eta_i.*(A./(GammaGth.*Lcav));
elseif isAVar && isEtaIVar
    eta_i=linspace(0.7, 0.9, 5)';
    A=linspace(0, 4, 32)';
    GammaGth=A./Lcav+alpha_i;

    for etaIndex=length(eta_i):-1:1
        %differential efficiency
        eta_d(:, etaIndex)=eta_i(etaIndex).*(A./(GammaGth.*Lcav));
    end
end



%%

figInterp=figure();
if isEtaIVar && ~isAVar
    plot(eta_d, eta_i, 'linewidth', 1.5)
    grid on
    xlabel('$\eta_d$','fontsize',30,'interpreter','latex')
    ylabel('$\eta_i$','fontsize',30,'interpreter','latex')
elseif isAVar && ~isEtaIVar
    plot(A, eta_d, 'linewidth', 1.5)
    grid on
    xlabel('A','fontsize',30,'interpreter','latex')
    ylabel('$\eta_d$','fontsize',30,'interpreter','latex')

elseif isAVar && isEtaIVar

    for etaIndex=length(eta_i):-1:1
        pl(etaIndex)=plot(A, eta_d(:, etaIndex), 'linewidth', 1.5, ...
            'DisplayName', ['$\eta_i$ = ', num2str(eta_i(etaIndex))]);
        hold on
        
    end

    grid on
    yline(0.2915, 'linestyle','--', 'linewidth', 1.5, ...
        'Label', ['$\eta_{d_{exp}} = $ ' , num2str(0.2915)], ...
        'fontsize',30, 'interpreter','latex')
    xlabel('A','fontsize',30,'interpreter','latex')
    ylabel('$\eta_d$','fontsize',30,'interpreter','latex')

    %legend settings
    legend(pl, 'location', 'southeast','NumColumns', 3,'fontsize',25,'interpreter','latex')
    box off
    legend boxoff
end

set(gca,'FontSize', 30,'TickLabelInterpreter','latex', 'TickDir', 'in') %Ticks

%figure settings
ScreenSize=get(groot,'ScreenSize');
ScreenSize=ScreenSize(3:4);
Height=600;
Width=1200;
set(gcf, 'color', 'w') ;
set(figInterp,'GraphicsSmoothing','on')
figInterp.Position=[(ScreenSize(1)-Width)/2 (ScreenSize(2)-Height)/2 Width Height];
