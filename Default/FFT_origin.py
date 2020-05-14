import numpy as np

import matplotlib.pyplot as plt

import csv

with open('outp.csv') as fp:
    original_data = list(csv.reader(fp))

M =len(original_data)
print(M)
time = np.zeros(M-1)
x_detec = np.zeros(M-1)
x_source = np.zeros(M-1)
for i in range(0,M-1):
    time[i] = original_data[i+1][1]
    x_source[i] = original_data[i+1][2]
    x_detec[i] = original_data[i+1][3]
print(time)
print(x_source)
print(x_detec)

#%matplotlib inline
plt.plot(time,x_detec)
plt.show()

plt.plot(time,x_source)
plt.show()

#%matplotlib inline
y = x_detec - x_detec[800]
print(y)
plt.plot(time,y)
plt.show()

c = 0
for i in range(M-1):
    if (i >=1000) and (y[i] <=0):
        u = y[i:i+600]
        #print(u)
        print(len(u))
        c = i + 1
        break
print(c)
#decide propagation time from the time when location becomes minus fitst,

import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack

xdata = time[c:c+600]
ydata = u 
print(len(xdata),len(ydata))
#plt.plot(xdata,ydata)
#plt.show()

time_step = xdata[1]-xdata[0]

#FFT
sample_freq = fftpack.fftfreq(ydata[:].size, d=time_step)
y_fft = fftpack.fft(ydata[:])
pidxs = np.where(sample_freq > 0)
freqs, power = sample_freq[pidxs], np.abs(y_fft)[pidxs]
freq = freqs[power.argmax()]

#PLot
plt.figure(figsize=(8,10))
plt.subplot(211)
plt.plot(xdata,ydata,'b-', linewidth=1)
plt.xlabel('Time')
plt.ylabel('Ydata')
plt.grid(True)

plt.subplot(212)
#plt.semilogx(freqs, power,'b.-',lw=1)
plt.loglog(freqs, power,'b.-',lw=1)
plt.xlabel('Frequency')
plt.ylabel('Power')
plt.grid(True)

plt.show()

print(power[2])
print(power[5])
print(power[8])

def calculate_beta(a1,a2,c):
    A1 = 0.05 * 2    #a1*2/600
    A2 = a2*2/600
    #v = 55.36231
    de =x_detec[800] -x_source[800]
    dti = (c-100)/100
    v = de/dti
    print(v)
    #pi = 3.141592
    beta = 8*A2*2*2*v*v/de/A1/A1/np.pi/np.pi/2/2
    return beta
#amplitude * 600

print(calculate_beta(power[2],power[5],c))
#error to the result of EXCEL is under 10^-5

x_bar = np.array([0.5,1,1.5,2,2.5,3])
y_bar = np.array([power[2],power[5],power[8],power[11],power[14],power[17]])
print(x_bar)
print(y_bar)
plt.bar(x_bar, y_bar,width=0.2)
plt.show()
#plt.plot()


















