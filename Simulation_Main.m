data = load('experiment_1_fan_on.mat');
data = data.DATAr;

fc = 76e9;
c = 3e8;
lambda = c/fc;
tm = 256e-6;
bw = 1e9;
sweep_slope = bw/tm;
fs = 2e6;

[r,a,z] = RAngCalc(data);


figure;
contourf(Ang,r,abs(t),20);
ylim([0.3 5])
caxis([0 2e4])
zlim([0 5e4])


function [Range, Angle, ReflecMap] = RAngCalc(data)
    tm = 256e-6;
    bw = 1e9;
    sweep_slope = bw/tm;
    fs = 2e6;
    
    ReflecMap = fftshift(fft2(data),2);
  
    Angle = acosd([-1 : 2/(size(data,2)-1) : 1]) - 90;
    Freq = [0 : fs/(size(data,1)-1) : fs];
    Range = (Freq*3e8) ./ (2*sweep_slope);

end





