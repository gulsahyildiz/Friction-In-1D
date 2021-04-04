import math
import matplotlib.pyplot as plt

print("------------------------------------------------------------------------------------------------------------------------------")
v0 = int(input("Initial Velocity (m/sec) : "))
h0 = int(input("Height (m) : "))

# Variables
p = 1.225 
A = 1 
Cd = 0.01 
k = (1/2)*p*A*Cd

m = 1  
g = 9.81

# Terminal Velocity
vt = math.sqrt(m*g/k)  

# Flight Times
t = 0
time_step = 0.0001
tmax = (vt/g)*math.atan(v0/vt)
ft = (vt/g)*math.acosh(math.exp(h0*g/(vt**2)))
ft2 = (vt/g)*(math.acosh(math.exp(g*h0/(vt**2) + math.log(abs(math.cosh(math.atan(v0/vt)))))) - math.atan(v0/vt))

# Decimal Range Function
def decimal_range(i, f, inc):
    while i < f:
        yield i
        i += inc

# Velocity & Height Calculations
time = []
velocity = []
height = []

if v0 == 0:
    for t in decimal_range(0, ft+time_step, time_step):
        t = t
        v = vt*math.tanh(t*(-g/vt))
        h = h0 - ((vt**2)/g)*math.log(math.cosh(t*-g/vt))
        if h >= 0:
                time.append(t)
                velocity.append(v)
                height.append(h)
else:
    if v0 < 0:
        for t in decimal_range(0, ft2+time_step, time_step):
            t = t
            v = vt*math.tanh(t*(-g/vt) + math.atanh(v0/vt))
            h = h0 + (vt**2/-g)*(math.log(abs(math.cosh(math.atan(v0/vt) - g*t/vt))) - math.log(abs(math.cosh(math.atan(v0/vt)))))
            if h >= 0:
                time.append(t)
                velocity.append(v)
                height.append(h)
    elif v0 > 0:
        for t in decimal_range(0, tmax, time_step):
            t = t
            v = vt*math.tan(math.atan(v0/vt) - g*t/vt)
            h = h0 - (m/k)*math.log((math.cos(math.atan(v0/vt)))/(math.cos(math.atan(v0/vt) - g*t/vt)))
            hmax = h
            ft0 = (vt/g)*math.acosh(math.exp(hmax*g/(vt**2)))
            time.append(t)
            velocity.append(v)
            height.append(h)
        for t in decimal_range(tmax, tmax+ft0+time_step, time_step):
            t = t
            v = vt*math.tanh((t - tmax)*(-g/vt))
            h = hmax - ((vt**2)/g)*math.log(math.cosh((t - tmax)*-g/vt))
            if h >= 0:
                time.append(t)
                velocity.append(v)
                height.append(h)

# Creating the S21.dat File
file = open("S21.dat", "w")
for n in range(len(time)):
  file.write("{0:0.4f} {1:30.4f} {2:60.4f}\n".format(time[n] , height[n], velocity[n]))
file.close()

# Reading the File
file = open("S21.dat", "r")
f = file.read()
print(f)

# Plotting the data
plt.plot(time, height, "r", label = "Height [m]") 
plt.plot(time, velocity, "b", label = "Velocity [m/s]")
plt.xlabel('Time')
plt.title('Height & Velocity vs Time')
plt.legend(loc = 1)
plt.minorticks_on()
plt.grid(b = True, which = 'both')
plt.show()

