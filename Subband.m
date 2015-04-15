%Prepare the Workspace
close all;
clear all;
clc;

%Load in Signal Samples
load money.dat;
x = subband';

%Display Signal
figure; plot(x, 'k');
title('Original Input Signal')
xlabel('Sample Number (n)');
ylabel('x(n)');
axis([0 length(x) -17000 22000]);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

%Display Impulse Response of H1
h1 = [0.0089874495 -0.0681367500 0.0669402090 0.4906729800 0.4906729800 0.0669402090 -0.0681367500 0.0089874495];

x1 = filter(h1, 1, x);
figure; plot(x1, 'r');
title('Low Pass Signal')
xlabel('Sample Number (n)');
ylabel('x1(n)');
axis([0 length(x) -17000 22000]);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

%Display Impulse Response of H2
h2 = h1; for i=2:2:length(h1) h2(i) = h2(i)*-1; end
x2 = filter(h2, 1, x);
figure; plot(x2, 'b');
title('High Pass Signal')
xlabel('Sample Number (n)');
ylabel('x2(n)');
axis([0 length(x) -17000 22000]);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

%Display Low & High Pass Filters with Combined Gain
figure; plot(linspace(0,1,1024),(abs(fft(h1,1024))),'k');
hold on; plot(linspace(0,1,1024),(abs(fft(h2,1024))),'r');
hold on; plot(linspace(0,1,1024),(abs(fft(h1+h2,1024))),'b--');
title('Magnitude Spectrum of Filters H1, H2 & H1+H2')
xlabel('Normalized Frequency (f)');
axis([0 0.5 0 1.1]);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

%Display the Combined Gain Error From Unity
figure; plot(linspace(0,1,1024), 20*log10(abs(abs(fft(h1,1024)+fft(h2,1024)))) ,'b');
title('Magnitude of Error dB of H1+H2')
xlabel('Normalized Frequency (f)');
ylabel('Error(f)');
axis([0 0.5 -0.03 0.03]);
grid on;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

k1 = h1.*2;   %Create Low-Pass Reconstruction Filter
k2 = h2.*-2;   %Create High-Pass Reconstruction Filter

%Down-Sample 2
x3 = downsample(x1, 2, 0);
x4 = downsample(x2, 2, 0);

%Cycle Thru Number of Bits From 1-to-7
for b=1:7;
    xlq=zeros(1,2000);
    x2q=zeros(1,2000);

    Qupper=max(x3)+10^-9;
    Qlower=min(x3);
    %Calculate Quanter Level
    delta = (Qupper-Qlower)/2^b;
    %Perform Quantization
    for n=1:length(x3)
       for i=Qlower:delta:Qupper
          if ((x3(n)>=i) && (x3(n)<=i+delta))  %Does point lie between?
             x1q(n) = i + delta / 2;    %Yes - Shift Value
          end
       end
    end
    
    Qupper=max(x4)+10^-9;
    Qlower=min(x4);
    %Calculate Quanter Level
    delta = (Qupper-Qlower)/2^(8-b);
    %Perform Quantization
    for n=1:length(x4)
       for i=Qlower:delta:Qupper
          if ((x4(n)>=i) && (x4(n)<=i+delta))  %Does point lie between?
             x2q(n) = i + delta / 2;    %Yes - Shift Value
          end
       end
    end
    
    %Up-Sample 2
    x3u = upsample(x1q, 2);
    x4u = upsample(x2q, 2);

    %Reconstruction Filters
    y1 = filter(k1, 1, x3u);
    y2 = filter(k2, 1, x4u);
    
    %Add the two reconstructed signals together
    y = y1 + y2;
    
    %Corrects the delay of e^-j*w*(M-1) introduced by the QMFs, this makes
    %for a realistic comparision between the input and output signal
    y = y([length(h1):end 1:length(h1)-1]);

    %If bits in low pass signal quantizer are equal to 5...
    if (b==5)
        figure; plot(x, 'r');
        hold on; plot(y, 'b');
        legend('Original Signal', 'Reconstructed Signal');
        xlabel('Sample Number (n)');
        ylabel('xhat(n)');
        axis([0 length(x) -17000 22000]);
        title('Original Signal vs. Reconstructed Signal');
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    end

    %SNR
    SNR(b) = 10 * log10( sum((x).^2) ./ sum((x-y).^2) );

end

% Compare SNR with PCM
figure; plot(SNR,'r'); hold on; grid on;
xlabel('Number of Bits in Low-Pass Filter'); ylabel('SNR');
title('SNR for Each Number of Bits vs. 4 bit PCM');
%Compute SNR for PCM
Qupper=max(x)+10^-9;
Qlower=min(x);
delta = (Qupper-Qlower)/2^4;   %4-Bit PCM System
%Perform Quantization
for n=1:length(x)
   for i=Qlower:delta:Qupper
      if ((x(n)>=i) && (x(n)<=i+delta))  %Does point lie between?
         xq(n) = i + delta / 2;    %Yes - Shift Value
      end
   end
end
eq = x-xq;
SNRpcm = 10*log10(sum(x.*x)/sum(eq.*eq));
line([1 7], [SNRpcm SNRpcm]);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
legend('SNR for Sub-Band Coding', 'SNR For PCM');
