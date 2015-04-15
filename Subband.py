import matplotlib.pyplot as plt
from scipy import signal
import scipy.fftpack as fp
import numpy as np

# Load in the signal
x = np.loadtxt('money.dat')

# Create QMF's
h1 = np.array([0.0089874495, -0.0681367500, 0.0669402090, 0.4906729800, 0.4906729800, 0.0669402090, -0.0681367500, 0.0089874495]);
x1 = signal.lfilter(h1, 1, x);
h2 = np.copy(h1)
for i in range(1, h2.__len__(), 2):
  h2[i] = h2[i] * -1;
x2 = signal.lfilter(h2, 1, x);
k1 = np.copy(h1*2);   #Create Low-Pass Reconstruction Filter
k2 = np.copy(h2*-2);   #Create High-Pass Reconstruction Filter

# Display the original signal
plt.figure()
plt.title('Original Input Signal')
plt.plot(x, 'k')
plt.axis([0, x.__len__(), -17000, 22000])
plt.xlabel('Sample Number (n)');
plt.ylabel('x(n)');

# Display Impulse Response of H1
plt.figure()
plt.plot(x1, 'r');
plt.axis([0, x.__len__(), -17000, 22000])
plt.title('Low Pass Signal')
plt.xlabel('Sample Number (n)');
plt.ylabel('x1(n)');

# Display Impulse Response of H2
plt.figure()
plt.plot(x2, 'b');
plt.title('High Pass Signal')
plt.xlabel('Sample Number (n)');
plt.ylabel('x2(n)');
plt.axis([0, x.__len__(), -17000, 22000])

# Display Magnitude Spectrum of Filters
plt.figure()
plt.hold(True)
plt.plot(np.arange(0., 1., 1/1024), np.abs(fp.fft(h1, 1024)), 'k');
plt.plot(np.arange(0., 1., 1/1024), np.abs(fp.fft(h2, 1024)), 'r');
plt.plot(np.arange(0., 1., 1/1024), np.abs(fp.fft(h1+h2, 1024)), 'b--');
plt.title('Magnitude Spectrum of Filters H1, H2 & H1+H2')
plt.xlabel('Normalized Frequency (f)');
plt.axis([0, 0.5, 0, 1.1]);

# Error between QMF filters
plt.figure()
plt.plot(np.arange(0., 1., 1/1024), 20*np.log10(np.abs(np.abs(fp.fft(h1,1024)+fp.fft(h2,1024)))) ,'b');
plt.title('Magnitude of Error dB of H1+H2')
plt.xlabel('Normalized Frequency (f)');
plt.ylabel('Error(f)');
plt.axis([0, 0.5, -0.03, 0.03]);

# Down-Sample 2 Matlab: downsample(x1, 2, 0)
x3 = x1[0::2]  # Format is [offset::downsample_factor]
x4 = x2[0::2]

# Cycle thru number of bits from 1-to-7
SNR = np.empty(8)
for b in np.arange(1,8):

  # Quantization 
  x1q = np.empty(x3.__len__())
  x2q = np.empty(x4.__len__())
  q1 = ((np.max(x3)+10**-9)-np.min(x3))/(2**b)
  q2 = ((np.max(x4)+10**-9)-np.min(x4))/(2**(8-b))
  x1q = np.min(x3) + np.floor((x3 - np.min(x3))/q1)*q1 + (q1/2)
  x2q = np.min(x4) + np.floor((x4 - np.min(x4))/q2)*q2 + (q2/2)

  # Plot the quantized signals
  #plt.figure()
  #plt.plot(x1q)
  #plt.figure()
  #plt.plot(x2q)

  # Upsample by 2
  x3u = np.append(np.insert(x1q, slice(1, x1q.__len__()), 0), [0])
  x4u = np.append(np.insert(x2q, slice(1, x2q.__len__()), 0), [0])

  # Filter
  y1 = signal.lfilter(k1, 1, x3u);
  y2 = signal.lfilter(k2, 1, x4u);

  # Add the two reconstructed signals together
  y = y1 + y2;

  # Corrects the delay of e^-j*w*(M-1) introduced by the QMFs, this makes
  # for a realistic comparision between the input and output signal
  y = np.append(y[h1.__len__()-1:y.__len__()],y[0:h1.__len__()-1])

  # When bit rate allocation is optimal, plot the reconstructed signal.
  if (b==5):
    plt.figure()
    plt.plot(x, "r")
    plt.hold(True)
    plt.plot(y, "b")
    plt.title('Reconstructed signal following Subband processing');
    plt.ylabel('y(n)');
    plt.xlabel('Sample Number');  

  # Calculate the SNR for each quantization bit rate
  SNR[b] = 10 * np.log10( np.sum((x)**2) / np.sum((x-y)**2) );

# Compare SNR with PCM
xq = np.empty(x.__len__())
q3 = ((np.max(x)+10**-9)-np.min(x))/(2**4)
xq = np.min(x) + np.floor((x - np.min(x))/q3)*q3 + (q3/2)
eq = x-xq
SNRpcm = [10 * np.log10( np.sum(x*x) / np.sum(eq*eq) )]*8;

# Display SNR comparison of subband vs PCM
plt.figure()
plt.plot(SNR, "b")
plt.hold(True)
plt.plot(SNRpcm, "r")
plt.axis([1, 7, -10, 20])
plt.title('Comparison of SNR of Subband vs PCM');
plt.ylabel('SNR (dB)');
plt.xlabel('Number of Bits');

# Show all the plots
plt.show()
