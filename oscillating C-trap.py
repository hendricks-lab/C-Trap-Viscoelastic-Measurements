import bluelake as bl
import numpy as np

freq = [0.097, 0.17, 0.27, 0.43, 0.57, 0.91, 1.23, 2.11, 3.03, 5.07, 8.77, 19.9, 36.7, 83, 170]
#freq =  [0.097,  0.17, 0.43,   1.37,    2.37,    5,      9.1,     19,     37,     87,     153]
#amp = [0.95,   0.9,    0.8,    0.4,     0.4,    0.29,   0.28,   0.25,   0.22,     0.2,    0.12]
amp = [0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.45, 0.4, 0.35, 0.29, 0.28, 0.25, 0.22, 0.2, 0.17]
amp = np.multiply(amp,10)
pauseTime = [30, 20, 10, 8, 7, 6, 5, 4, 3, 3, 2, 1, 1, 1, 1]

print(pauseTime)
x, y = bl.mirror1.position
for i in range(0,len(freq)):
    try:
        bl.mirror1.move_to(x=x, y=y, speed=10)
        print(bl.mirror1.position)
        bl.mirror1.start_oscillation(axis='X',amplitude = amp[i],frequency = freq[i])
        bl.pause(pauseTime[i])

    finally:
        bl.mirror1.stop()


limits = bl.mirror1.motion_limits

# Print coordinates
print(limits.lower.x)
print(limits.upper.x)