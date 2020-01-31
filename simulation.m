%%
% The following table summarizes the radar parameters.
% 
%  System parameters            Value
%  ----------------------------------
%  Operating frequency (GHz)    77
%  Sweep time (microseconds)    60
%  Sweep bandwidth (GHz)        2
%  Sample rate (GHz)            2

fc = 77e9;
c = 3e8;
lambda = c/fc;
tm = 60e-6;
bw = 2e9;
sweep_slope = bw/tm;
fs = 2*bw;   %sample rate

%% Create Waveform
waveform = phased.FMCWWaveform('SweepTime',tm,'SweepBandwidth',bw,...
    'SampleRate',fs);


%% create patch antenna array
%Define Array
Ant =  phased.CosineAntennaElement;
Ant.FrequencyRange = [fc fc+bw];

AntArray = phased.ULA;
AntArray.Element = Ant;
AntArray.NumElements = 8;
AntArray.ElementSpacing = 0.5*lambda;

%% Steering Pattern By 10 degrees
%steervec = phased.SteeringVector('SensorArray',AntArray);
%sv = steervec(fc,[10;0]);
subplot(131);
pattern(AntArray,fc)
%title('Without steering')
%subplot(122)
%pattern(AntArray,fc,'Weights',sv)

%% Transmitter and Receiver
tx_ppower = db2pow(10)*1e3;            % in watts
tx_gain = 17;                          % in dB

rx_gain = 15;                          % in dB
rx_nf = 4.5;                           % in dB

transmitter = phased.Transmitter('PeakPower',tx_ppower,'Gain',tx_gain);
receiver = phased.ReceiverPreamp('Gain',rx_gain,'NoiseFigure',rx_nf,...
    'SampleRate',fs);

%% Radiator
radiator = phased.Radiator('Sensor',AntArray,'OperatingFrequency',fc);

%% Collector
collector = phased.Collector('Sensor',AntArray,'OperatingFrequency',fc);

%% Target and Stuff
tgt_pos = [5;0;0];
tgt_vel = [0;0;0];
tar_rcs = 0.01;

target = phased.RadarTarget('MeanRCS',tar_rcs,'PropagationSpeed',c,...
    'OperatingFrequency',fc);
tarmotion = phased.Platform('InitialPosition',tgt_pos,...
    'Velocity',tgt_vel);

radar_pos = [0;0;0];
radar_vel = [0;0;0];
radarmotion = phased.Platform('InitialPosition',radar_pos,...
    'Velocity',radar_vel);


%% Propagation

channel = phased.FreeSpace(...
    'PropagationSpeed',physconst('LightSpeed'),...
    'OperatingFrequency',fc,'TwoWayPropagation',false,...
    'SampleRate',fs);


%% Loop

Nsweep = 64;
xr = complex(zeros(waveform.SampleRate*waveform.SweepTime,Nsweep));

for m = 1:Nsweep
    
    % Transmit FMCW waveform
    sig = waveform();
    txsig = transmitter(sig);
    
    txsig = radiator(txsig,0);

    % Propagate the signal and reflect off the target
    txsig = channel(txsig,radar_pos,tgt_pos,radar_vel,tgt_vel);
    txsig = target(txsig);
    txsig = channel(txsig,tgt_pos,radar_pos,radar_vel,tgt_vel);
    % Dechirp the received radar return
    txsig = collector(txsig,0);
    txsig = receiver(txsig);
    
    %dechirpsig = dechirp(txsig,sig);

    % Visualize the spectrum
    %specanalyzer([txsig dechirpsig]);

    %xr(:,m) = dechirpsig;
end


subplot(132); plot(0:1/fs:tm-1/fs,real(txsig(:,1)));
subplot(133); plot(0:1/fs:tm-1/fs,real(txsig(:,8)));
xlabel('Time (s)'); ylabel('Amplitude (v)');
title('FMCW signal'); axis tight;
%subplot(212); spectrogram(sig,32,16,32,fs,'yaxis');
title('FMCW signal spectrogram');







