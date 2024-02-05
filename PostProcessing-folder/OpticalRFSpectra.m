function OpticalRFSpectra(Field, TStart, TEnd)
% TStart=2;
% TEnd=Inf; %why?? 
Subsample=1;
color='b';
V=2240;
L=2000; 

%% Figure settings
ScreenSize=get(groot,'ScreenSize');
ScreenSize=ScreenSize(3:4);

Height=600;
Width=1300;

figSpec=figure();
                               
%% Generation of the plot of Optical spectrum
valid=(Field.time>=TStart & Field.time<=TEnd).*(rem(1:length(Field.time), Subsample)==0)==1;
t=Field.time(valid);
y=Field.EL(valid);
Y=fftshift(fft(y))/length(y);
dt=t(2)-t(1);
f=linspace(-0.5/dt,+0.5/dt,length(t));
P=abs(Y).^2;
subplot(1,2,1)
plot(f,10*log10(P),'color', [0.39,0.83,0.07]);
xlim([-600 600])
ylim([-110 10])
hold on;
%
xlabel('Frequency (GHz)','fontsize',30,'interpreter','latex')
ylabel('Optical Spectrum (dB)','fontsize',30,'interpreter','latex')

set(gca,'xminortick', 'on', ...
    'yminortick', 'on','FontSize', 30,...
        'TickLabelInterpreter','latex', 'TickDir', 'in') %Ticks

%% Generation of the plot of the Power spectrum 
subplot(1,2,2)
valid=(Field.time>=TStart & Field.time<=TEnd).*(rem(1:length(Field.time), Subsample)==0)==1;
t=Field.time(valid);
y=abs(Field.EL(valid)).^2;
Y=fftshift(fft(y))/length(y);
dt=t(2)-t(1);
f=linspace(-0.5/dt,+0.5/dt,length(t));
Pf=abs(Y).^2;
plot(f,10*log10(Pf),'color', [0.12,0.01,0.32] );
xlim([0 200])
ylim([-170 0])
xlabel('Frequency (GHz)','fontsize',30,'interpreter','latex')
ylabel('Power Spectrum (dB)','fontsize',30,'interpreter','latex')

set(gca,'xminortick', 'on', ...
    'yminortick', 'on','FontSize', 30,...
        'TickLabelInterpreter','latex', 'TickDir', 'in') %Ticks

set(figSpec,'GraphicsSmoothing','on')
set(gcf, 'color', 'w') ;
figSpec.Position=[(ScreenSize(1)-Width)/2 (ScreenSize(2)-Height)/2 Width Height];