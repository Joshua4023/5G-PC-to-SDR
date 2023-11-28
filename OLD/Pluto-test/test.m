
% Number of points in the waveform
points = 10000;

% Determines the frequency offset from the carrier
cycles = 101;
phaseInc = 2*pi*cycles/points;
phase = phaseInc * (0:points-1);

% Create an IQ waveform
Iwave = cos(phase);
Qwave = sin(phase);
IQData = Iwave+1i*Qwave;
IQData = IQData(:)';
IQData = rot90(IQData);





%% Visualize

time = 10*1e-3;
length = size(IQData);

x = 0:65e-9:time;
x = rot90(x);

x=x(1:length, 1);





plot(x, real(IQData),x, imag(IQData));
legend("Real", "Imag")
xlabel('Time(S)')


%% Transmit waveform over the air
plutoTx = sdrtx('Pluto', RadioID='usb:0');
plutoTx.CenterFrequency = 75e6;
plutoTx.Gain = -80;
plutoTx.BasebandSampleRate = 61440000;
plutoTx.ShowAdvancedProperties = true;
plutoTx.FrequencyCorrection = 0;

% Transmit waveform (for 10 sec):
transmitRepeat(plutoTx, IQData);