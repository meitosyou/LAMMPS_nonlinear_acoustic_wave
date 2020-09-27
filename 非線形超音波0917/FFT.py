import numpy as np
# グラフを描画するためのもの
import matplotlib.pyplot as plt
import csv

with open('outp.csv') as fp:
    original_data = list(csv.reader(fp))

#自分で設定した系のものと合わせる。
A1 = .10
period = 2
timestep = 0.001
T = int(period/timestep*10)
T_wave = int(period/timestep*4)

M =len(original_data)
print(M)
time = np.zeros(M-1)
x_detec = np.zeros(M-1)
x_source = np.zeros(M-1)

for i in range(0,M-1):
    time[i] = original_data[i+1][0]
    x_source[i] = original_data[i+1][1]
    x_detec[i] = original_data[i+1][2]
print(time)
print(x_source)
print(x_detec)
"""
%matplotlib inline
plt.plot(time,x_detec)
plt.show()

plt.plot(time,x_source)
plt.show()

%matplotlib inline
y = x_detec - x_detec[T_wave]
print(y)
plt.plot(time,y)
plt.show()
"""

y = x_detec - x_detec[T_wave]
print(y)

#到達時間判定法　detectorの位置が上がってさがって、初めてもとの位置を下回ったときが半周期に当たるはずだから、その時間から半周気分を引いた時間を到達時間としている。この時間判定法は時間ステップの粗さなどから少々課題がある。
c = 0
for i in range(M-1):
    if (i >=int(period/timestep*5)) and (y[i] <=0):
        u = y[i:i+int(period/timestep*3)]
        #print(u)
        print(len(u))
        c = i + 1
        break
print(c)

#FFTをする部分
import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack

xdata = time[c:c+int(period/timestep*3)]
ydata = u 
print(len(xdata),len(ydata))
#plt.plot(xdata,ydata)
#plt.show()

#time_step = xdata[1]-xdata[0]
time_step = timestep
#FFT
sample_freq = fftpack.fftfreq(ydata[:].size, d=time_step)
y_fft = fftpack.fft(ydata[:])
pidxs = np.where(sample_freq > 0)
freqs, power = sample_freq[pidxs], np.abs(y_fft)[pidxs]
freq = freqs[power.argmax()]

#plot
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


#基本波、高調波の振幅を見てみる
print(power[2])
print(power[5])
print(power[8])

#非線形パラメータbetaを算出する。
def calculate_beta(a1,a2,c):
    print(A1)
    A2 = a2*2/int(period/timestep*3)
    #print(int(period/timestep*3))
    print(A2)
    #v = 55.36231
    de =x_detec[T_wave] -x_source[T_wave]
    print(de)
    dti = (c-period/timestep/2)*timestep
    v = de/dti
    print(v)
    print("lambda")
    print(v*period)
    #pi = 3.141592
    beta = 8*A2*period*period*v*v/de/A1/A1/np.pi/np.pi/2/2
    return beta

#なんか振幅へん。なんで30?サンプル数の６００で割ればいい？
#棒グラフ書いてまとめようか
print("CALCULATE_BETA")
print(calculate_beta(power[2],power[5],c))
print("CALCULATE_COMPLETED")
#EXCELの結果と誤差が10^-5以下

x_bar = np.array([0.5,1,1.5,2,2.5,3])
y_bar = np.array([power[2],power[5],power[8],power[11],power[14],power[17]])
print(x_bar)
print(y_bar)

print(c)
print(power[2]/int(period/timestep*3)*2)
print(power[5]/int(period/timestep*3)*2)

plt.bar(x_bar, y_bar,width=0.2)
plt.show()


####file output
#file = "result.txt"
#fileobj = open(file, "w", encoding = "utf_8")
#fileobj.write("こんにちは\n")
#fileobj.close()
#with open('result.txt', 'w') as f:
#  print(A1,power[5]/int(period/timestep*3)*2 , file=f)

