%losses=Field.Losses;
function post_proc_simple_NEW(Field, Carriers, TStart, TEnd)
% TStart=2;
% TEnd=Inf; %why?? 
Subsample=1;
color='b';
V=2400;
L=2000; 

%% Figure settings
ScreenSize=get(groot,'ScreenSize');
ScreenSize=ScreenSize(3:4);

Height=600;
Width=1200;

%% Generation of the plot of Power vs time
valid=(Field.time>=TStart & Field.time<=TEnd).*(rem(1:length(Field.time), Subsample)==0)==1;
t=Field.time(valid);
y=Field.EL(valid);
%what are these numbers
%PLtime=0.5*0.127*2.997925e5*3.6*8.854e-18/(1.6e-19)*(1-0.3)*14*(2/25)*abs(y).^2*1.6*10^(-7);%.*1.6e-19*1/ResTD.Sim.dt/1e-9.*1e3;
PLtime=(0.5*2.997925e5*3.3*8.854e-18/(1.6e-19)*...
            (abs(y).^2)*(1-0.3)*(1/0.127)*V/L)*(1.6e-8);

figure ('Position',[(ScreenSize(1)-Width)/2 (ScreenSize(2)-Height)/2 Width Height])

%title(titlestr, 'fontsize',30,'interpreter','latex')
subplot(2,3,1)
plot(t,PLtime,'color',color); 
hold on; 
grid on;
xlabel('Time [ns]');
ylabel('Output power [mW]');
                               
%% Generation of the plot of Optical spectrum
valid=(Field.time>=TStart & Field.time<=TEnd).*(rem(1:length(Field.time), Subsample)==0)==1;
t=Field.time(valid);
y=Field.EL(valid);
Y=fftshift(fft(y))/length(y);
dt=t(2)-t(1);
f=linspace(-0.5/dt,+0.5/dt,length(t));
P=abs(Y).^2;
subplot(2,3,2)
plot(f,P,'color','r');
hold on;
%
xlabel('Frequency [GHz]')
ylabel(' Optical Spectrum [dB]')

%% Generation of the plot of the Power spectrum 
subplot(2,3,3)
hold on; 
valid=(Field.time>=TStart & Field.time<=TEnd).*(rem(1:length(Field.time), Subsample)==0)==1;
t=Field.time(valid);
y=abs(Field.EL(valid)).^2;

Y=fftshift(fft(y))/length(y);
dt=t(2)-t(1);
f=linspace(-0.5/dt,+0.5/dt,length(t));
Pf=abs(Y).^2;
plot(f,10*log10(Pf),'color','k');
xlim([0 100])
hold on;

xlabel('Frequency [GHz]')
ylabel(' Power Spectrum [dB]')

%% Generation of the plot of the Carrier dynamics
valid=(Field.time>=TStart & Field.time<=TEnd).*(rem(1:length(Field.time), Subsample)==0)==1;
t=Field.time(valid);

y0=Carriers.rhoGS0(:,valid)';
yp=Carriers.rhoGSp(:,valid)';

ypm=mean(abs(yp));
sprintf('mean %g yp',ypm)

subplot(2,3,4)
plot(t,y0); 
grid on;
xlabel('Time [ns]');
ylabel('Carriers order 0 (Z=L)');

subplot(2,3,5)
plot(t,2*real(yp)); 
grid on;
%title(sprintf('From %g to %g ns',t(1),t(end)));
xlabel('Time [ns]');
ylabel('Carriers grating (Z=L)');

%Plot of the instantaneous frequency
valid=(Field.time>=TStart & Field.time<=TEnd).*(rem(1:length(Field.time), Subsample)==0)==1;
t=Field.time(valid);
y=Field.EL(valid);
fieldphase=unwrap(angle(y));
inst=(1/(2*pi)).*diff(fieldphase)./(t(2)-t(1));
%lambdainst=((c/n_index)./inst)*1e-6;
subplot(2,3,6)
plot(t(1:length(inst)),inst); 
grid on;

xlabel('Time [ns]');
ylabel('Instantaneous Frequency [GHz]');

set(gcf, 'color', 'w');

%% Simultaneous plot of Power and Instantaneous frequency

figure
yyaxis left
plot(t,PLtime,'color',color); 
xlabel('Time [ns]');
ylabel('Output power [mW]');

yyaxis right
plot(t(1:length(inst)),inst); 
xlabel('Time [ns]');
ylabel('Instantaneous Frequency [GHz]')

set(gcf, 'color', 'w');






