fc = 77e9;
c = 3e8;
lambda = c/fc;
tm = 60e-6;
bw = 2e9;
sweep_slope = bw/tm;

v_max = 5;
fs = bw;

%%
% The following table summarizes the radar parameters.
% 
%  System parameters            Value
%  ----------------------------------
%  Operating frequency (GHz)    77
%  Sweep time (microseconds)    60
%  Sweep bandwidth (GHz)        2
%  Sample rate (GHz)            2

waveform = phased.FMCWWaveform('SweepTime',tm,'SweepBandwidth',bw,...
    'SampleRate',fs);

tar_dist = 20;
tar_speed = 0;
tar_rcs = 0.01;

target = phased.RadarTarget('MeanRCS',tar_rcs,'PropagationSpeed',c,...
    'OperatingFrequency',fc);
tarmotion = phased.Platform('InitialPosition',[tar_dist;0;0.5],...
    'Velocity',[tar_speed;0;0]);

radar_speed = 0;
radarmotion = phased.Platform('InitialPosition',[0;0;0.5],...
    'Velocity',[radar_speed;0;0]);

channel = phased.FreeSpace('PropagationSpeed',c,...
    'OperatingFrequency',fc,'SampleRate',fs,'TwoWayPropagation',true);

tx_ppower = db2pow(10)*1e3;            % in watts
tx_gain = 17;                          % in dB

rx_gain = 15;                          % in dB
rx_nf = 4.5;                           % in dB

transmitter = phased.Transmitter('PeakPower',tx_ppower,'Gain',tx_gain);
receiver = phased.ReceiverPreamp('Gain',rx_gain,'NoiseFigure',rx_nf,...
    'SampleRate',fs);


Nsweep = 64;
xr = complex(zeros(waveform.SampleRate*waveform.SweepTime,Nsweep));

for m = 1:Nsweep
    % Update radar and target positions
    [radar_pos,radar_vel] = radarmotion(waveform.SweepTime);
    [tgt_pos,tgt_vel] = tarmotion(waveform.SweepTime);

    % Transmit FMCW waveform
    sig = waveform();
    txsig = transmitter(sig);
    
    % Propagate the signal and reflect off the target
    txsig = channel(txsig,radar_pos,tgt_pos,radar_vel,tgt_vel);
    txsig = target(txsig);
    
    % Dechirp the received radar return
    txsig = receiver(txsig); 
  
    dechirpsig = dechirp(txsig,sig);
    
    % Visualize the spectrum
    specanalyzer([txsig dechirpsig]);
    
    xr(:,m) = dechirpsig;
end

sig = waveform();
subplot(211); plot(0:1/fs:tm-1/fs,real(txsig));
xlabel('Time (s)'); ylabel('Amplitude (v)');
title('FMCW signal'); axis tight;
subplot(212); spectrogram(sig,32,16,32,fs,'yaxis');
title('FMCW signal spectrogram');


rngdopresp = phased.RangeDopplerResponse('PropagationSpeed',c,...
    'DopplerOutput','Speed','OperatingFrequency',fc,'SampleRate',fs,...
    'RangeMethod','FFT','SweepSlope',sweep_slope,...
    'RangeFFTLengthSource','Property','RangeFFTLength',2048,...
    'DopplerFFTLengthSource','Property','DopplerFFTLength',256);

clf;
plotResponse(rngdopresp,xr);                     % Plot range Doppler map
axis([-v_max v_max 0 50])
clim = caxis;
