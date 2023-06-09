import bluelake as bl
import numpy as np
from bluelake import mirror1

freq = [1, 5, 10];
amp = [1,0.4,0.1];
amp = np.multiply(amp,7);

dif_amp = np.abs(np.diff(amp));
dif_amp = np.append(dif_amp,0);
print(dif_amp)

pauseTime = [10, 10, 10];
center_x, center_y = bl.mirror1.position

for i in range(0,len(freq)):
    try:

        bl.mirror1.start_oscillation(axis='X',amplitude = amp[i],frequency = freq[i])
        bl.pause(pauseTime[i])
        bl.mirror1.stop()
        bl.mirror1.move_to(x= center_x, y = center_y, speed =10)
        bl.mirror1.move_by(dx = dif_amp[i], dy = 0, speed =10)
        center_x, center_y = bl.mirror1.position

    finally:
        bl.mirror1.stop()




# Print coordinates
print('DONE')