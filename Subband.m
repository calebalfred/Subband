%% Prepare the Workspace
close all;
clear all;
clc;

% savePlots = true - To save the graphs to file, no viewing in matlab.
savePlots = false;

%% Load in Signal Samples
load money.dat;  % Headerless audio file saying "Money"
x = money';

%% Display the Original Signal
figure; plot(x, 'k');
title('Original Input Signal')
xlabel('Sample Number (n)');
ylabel('x(n)');
axis([0 length(x) -17000 22000]);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
if savePlots == true, saveas(gcf, '1 - Original Signal', 'png'), end

%% Display Impulse Response of H1
h1 = [0.0089874495 -0.0681367500 0.0669402090 0.4906729800 0.4906729800 0.0669402090 -0.0681367500 0.0089874495];

x1 = filter(h1, 1, x);
figure; plot(x1, 'r');
title('Low Pass Signal')
xlabel('Sample Number (n)');
ylabel('x1(n)');
axis([0 length(x) -17000 22000]);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
if savePlots == true, saveas(gcf, '2 - Low Passed Signal', 'png'), end

%% Display Impulse Response of H2
h2 = h1; for i=2:2:length(h1) h2(i) = h2(i)*-1; end
x2 = filter(h2, 1, x);
figure; plot(x2, 'b');
title('High Pass Signal')
xlabel('Sample Number (n)');
ylabel('x2(n)');
axis([0 length(x) -17000 22000]);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
if savePlots == true, saveas(gcf, '3 - High Passed Signal', 'png'), end

%% Display Low & High Pass Filters with Combined Gain
figure; plot(linspace(0,1,1024),(abs(fft(h1,1024))),'k');
hold on; plot(linspace(0,1,1024),(abs(fft(h2,1024))),'r');
hold on; plot(linspace(0,1,1024),(abs(fft(h1+h2,1024))),'b--');
title('Magnitude Spectrum of Filters H1, H2 & H1+H2')
xlabel('Normalized Frequency (f)');
axis([0 0.5 0 1.1]);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
if savePlots == true, saveas(gcf, '4 - H & L Filters Combined Gain', 'png'), end

%% Display the Combined Gain Error From Unity
figure; plot(linspace(0,1,1024), 20*log10(abs(abs(fft(h1,1024)+fft(h2,1024)))) ,'b');
title('Magnitude of Error dB of H1+H2')
xlabel('Normalized Frequency (f)');
ylabel('Error(f)');
axis([0 0.5 -0.03 0.03]);
grid on;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
if savePlots == true, saveas(gcf, '5 - H & L Filters Combined Error', 'png'), end

%% Create the reconstruction filters
k1 = h1.*2;   % Create Low-Pass Reconstruction Filter
k2 = h2.*-2;   % Create High-Pass Reconstruction Filter

%% Down-Sample 2
x3 = downsample(x1, 2, 0);
x4 = downsample(x2, 2, 0);

%% Initializing 
x1q = zeros(7, 2000);
x2q = zeros(7, 2000);
x3u = zeros(7, 4000);
x4u = zeros(7, 4000);
y1 = zeros(7, 4000);
y2 = zeros(7, 4000);
y = zeros(7, 4000);

%% Cycle Thru Number of Bits From 1-to-7
for b=1:7
    % Quantize
    q1(b) = ((max(x3)+10^-9)-min(x3))/(2^b);
    q2(b) = ((max(x4)+10^-9)-min(x4))/(2^(8-b));
    x1q(b, :) = min(x3) + floor((x3 - min(x3))/q1(b))*q1(b) + (q1(b)/2);
    x2q(b, :) = min(x4) + floor((x4 - min(x4))/q2(b))*q2(b) + (q2(b)/2);
    
    % Up-Sample 2
    x3u(b, :) = upsample(x1q(b, :), 2);
    x4u(b, :) = upsample(x2q(b, :), 2);

    % Reconstruction Filters
    y1(b, :) = filter(k1, 1, x3u(b, :));
    y2(b, :) = filter(k2, 1, x4u(b, :));
    
    % Add the two reconstructed signals together
    y(b, :) = y1(b, :) + y2(b, :);
    
    % Corrects the delay of e^-j*w*(M-1) introduced by the QMFs, this makes
    % for a realistic comparision between the input and output signal
    y(b, :) = y(b, [length(h1):end 1:length(h1)-1]);

    % SNR
    SNR(b) = 10 * log10( sum((x).^2) ./ sum((x-y(b, :)).^2) );
end

%% Display Quantization Plots
for b = 1:7
    figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    title('Quantization Bit Distribution - Low vs High ');
    subplot(2, 1, 1);
    plot(x1q(b, :), 'r');
    title( sprintf('Low Frequency - %d-bits', b) )
    subplot(2, 1, 2);
    plot(x2q(b, :), 'b');
    title( sprintf( 'High Frequency - %d-bits', (8-b) ) )
    if savePlots == true, saveas(gcf, sprintf( '%d - Quantization Distribution - L%d-H%d', 5+b, b, (8-b)), 'png'), end
end

%% Plot the reconstructed signal
figure; plot(x, 'r');
hold on; plot(y(5, :), 'b');
legend('Original Signal', 'Reconstructed Signal');
xlabel('Sample Number (n)');
ylabel('xhat(n)');
axis([0 length(x) -17000 22000]);
title('Original Signal vs. Reconstructed Signal');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
if savePlots == true, saveas(gcf, '13 - Reconstructed Signal', 'png'), end

%% Compare SNR with PCM
% Quantize the original signal using pcm quality
q3 = ((max(x)+10^-9)-min(x))/(2^4);
x3q = min(x) + floor((x - min(x))/q3)*q3 + (q3/2);
% Calculate the SNR of PCM quality signal
eq = x-x3q;
SNRpcm = 10*log10(sum(x.*x)/sum(eq.*eq));
% Plot the results
figure; plot(SNR,'r'); hold on; grid on;
line([1 7], [SNRpcm SNRpcm]);
legend('SNR for Sub-Band Coding', 'SNR For PCM');
xlabel('Number of Bits in Low-Pass Filter'); ylabel('SNR');
title('SNR for Each Number of Bits vs. 4 bit PCM');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
if savePlots == true, saveas(gcf, '14 - Bits 5L-3H SNR vs PCM', 'png'), end

%% Close down
if savePlots == true 
   close all;
end
clear all;
clc;
