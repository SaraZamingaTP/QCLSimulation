%losses=Field.Losses;
T0=2;
TEnd=Inf;
Subsample=1;
color='b';

%% Generation of the plot of Power vs time
valid=(Field.time>=T0 & Field.time<=TEnd).*(rem(1:length(Field.time), Subsample)==0)==1;
t=Field.time(valid);
y=Field.EL(valid);
PLtime=0.5*0.127*2.997925e5*3.6*8.854e-18/(1.6e-19)*(1-0.3)*14*(2/25)*abs(y).^2*1.6*10^(-7);%.*1.6e-19*1/ResTD.Sim.dt/1e-9.*1e3;
PLtimeAv=mean(PLtime(t>t(end)-1));
figure ('Position',[260 278 960 620])
subplot(2,3,1)
plot(t,PLtime,'color',color); 
hold on; 
grid on;
xlabel('Time [ns]');
ylabel('Output power [mW]');

hPlack=4.1357e-6;     %[ev ns]                                     %Plank constant; [eV ns]
PLfield=abs(y);
PLfieldAv=mean(PLfield(t>t(end)-1));


%% Generation of the plot of Optical spectrum
valid=(Field.time>=T0 & Field.time<=TEnd).*(rem(1:length(Field.time), Subsample)==0)==1;
t=Field.time(valid);
y=Field.EL(valid);
Y=fftshift(fft(y))/length(y);
dt=t(2)-t(1);
f=linspace(-0.5/dt,+0.5/dt,length(t));
lambda=1.24./(f*hPlack)+10; 
omega=(f*2*pi); %angular frequency in GHz; 
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
valid=(Field.time>=T0 & Field.time<=TEnd).*(rem(1:length(Field.time), Subsample)==0)==1;
t=Field.time(valid);
y=abs(Field.EL(valid)).^2;

Y=fftshift(fft(y))/length(y);
dt=t(2)-t(1);
f=linspace(-0.5/dt,+0.5/dt,length(t));
%lambda=1.24./(EnergyGapQDcentral+f*hPlack)*1000;
omega=(f*2*pi); %angular frequency in GHz; 
Pf=abs(Y).^2;
valid1=(omega>0.1);
Ypos=10*log10(Pf(valid1));
omegapos=omega(valid1);
[maxy,imax]=max(Ypos);
omegafmax=omegapos(imax);
plot(f,10*log10(Pf),'color','k');
xlim([0 100])
hold on;
%
xlabel('Frequency [GHz]')
ylabel(' Power Spectrum [dB]')

%% Generation of the plot of the Carrier dynamics
valid=(Field.time>=T0 & Field.time<=TEnd).*(rem(1:length(Field.time), Subsample)==0)==1;
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
title(sprintf('From %g to %g ns',t(1),t(end)));
xlabel('Time [ns]');
ylabel('Carriers grating (Z=L)');

%Plot of the instantaneous frequency
valid=(Field.time>=T0 & Field.time<=TEnd).*(rem(1:length(Field.time), Subsample)==0)==1;
t=Field.time(valid);
y=Field.EL(valid);
fieldphase=unwrap(angle(y));
inst=(1/(2*pi)).*diff(fieldphase)./(t(2)-t(1));
subplot(2,3,6)



plot(t(1:length(inst)),inst); 
grid on;

xlabel('Time [ns]');
ylabel('Instantaneous Frequency [GHz]');

%Simultaneous plot of Power and Instantaneous frequency

figure
yyaxis left
plot(t,PLtime,'color',color); 
xlabel('Time [ns]');
ylabel('Output power [mW]');

yyaxis right
plot(t(1:length(inst)),inst); 
xlabel('Time [ns]');
ylabel('Instantaneous Frequency [GHz]')


