data = load('experiment_1_fan_on.mat');
data = data.DATAr;

fc = 76e9;
c = 3e8;
lambda = c/fc;
tm = 256e-6;
bw = 1e9;
k = bw/tm;
fs = 2e6;   %sample rate after stretch proccessing
fss = 2e9;  %sample rate to generate waveforms
t = [0 : 1/fss : tm];

%% Geom Sim 2D
% assume Tx element is a origin but assuming target radiates isotropically
% atm so doesnt matter
RxArray = zeros(29,2);
RxSig = zeros((length(t)-1)/(fss/fs),29);
Reflector = [0,5];


for i = 1:29
    RxArray(i,:) = [i 0] * (lambda/2); %seperate 29 elements by distance lambda/2 
    tdelay = (sqrt((Reflector(1)-RxArray(i,1))^2 + (Reflector(2)-RxArray(i,1))^2)/c);
    s_delay = exp(2*pi*1j*(fc*(t-tdelay)+k*(t-tdelay).^2));
    mix = exp(2*pi*1j*(fc*t+k*t.^2)) .* conj(s_delay);
    RxSig(:,i) = mix(1:(fss/fs):(end-1));
end

[rs,as,Ran,Ang,zs] = RAngCalc(RxSig,k,fs);

figure;
subplot(131);
for i=1:29
    plot(rs,abs(Ran(:,i)));
    hold on
end

subplot(132);
for i=1:29
    plot(as,abs(Ang(i,:)));
    hold on
end
subplot(133);
contourf(as,rs,abs(zs))
ylim([0 10])
disp(2);



%% Real Data Bit
[r,a,Ra,An,z] = RAngCalc(data,k,fs);

figure;
subplot(131);
for i=1:29
    plot(r,abs(Ra(:,i)));
    hold on
end

subplot(132);
for i=1:29
    plot(a,abs(An(i,:)));
    hold on
end

subplot(133);
contourf(a,r,abs(z),20);
ylim([0.3 5])
caxis([0 2e4])
zlim([0 5e4])


%% Functions
function [Range, Angle, R, A, ReflecMap] = RAngCalc(data,sweep_slope,fs) 
 
    R = fft(data,[],1);
    A = fftshift(fft(data,[],2),2);
    ReflecMap = fftshift(fft(R,[],2),2);
    Angle = asind([-1 : 2/(size(data,2)-1) : 1]);
    Freq = [0 : fs/(size(data,1)-1) : fs];
    Range = (Freq*3e8) ./ (2*sweep_slope);

end





