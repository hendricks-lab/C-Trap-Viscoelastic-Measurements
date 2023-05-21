#!/usr/bin/env python
# coding: utf-8

# In[4]:


# Ben 2023 - CTRAP automated trap steering
# This code takes a variety of waveforms - currently sin, square, and triangle waveforms and breaks them down into discrete time steps.
# The CTrap is cappable of passing simple sinusoidal oscillations - but there are no internal commands to pass other types of trap displacements.
# This code takes whatever waveform and breaks it down into individual steps that the trap mirros can do using the bl.mirror1.move_by command.
# Set the desired displacement to True and run code.

import lumicks.pylake as lk
import numpy as np
import matplotlib.pyplot as plt
import time
#import bluelake as bl

# NOTE: ADD ANY DISPLACEMENT WAVEFORMS IN THIS INITIAL SECTION - MAKE SURE TO INCLUDE THE FOLLOWING LINES OF CODE AT THE END OF THE DISPLACEMENT DEFINITION:
    # function = [the y values that are used to plot the displacement]
    # num_points = len(function)
# THESE ARE USED WHEN CREATING THE INDIVIDUAL DISPLACEMENT STEPS. 

# Incremental Sinusoidal Waveform Generation
    # sampling_rate--> the number of point per second used to generate the square wave
if False:
    amplitude = 2
    frequency = 2
    phase = 0 #np.pi/4
    duration = 1
    sampling_rate = 1000
    time_array = np.arange(0, duration, 1/sampling_rate)

    sinusoidal_function = amplitude * np.sin(2*np.pi*frequency*time_array + phase) # the number of discrete point thaat make up this function is the time * the sampling rate. 

    plt.plot(time_array, sinusoidal_function) # the positions corresponding to the 
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title('Continuous Sinusoidal Function')
    plt.show()
    
    function = sinusoidal_function;
    num_points = len(sinusoidal_function); # essentially equivalent to the sampling rate.
  
# Oscillation Sinusoidal Waveform Generation
    # just uses the regular bluelake commands to pass oscillations - not incremental.
if False:  
    A = 0.5; # (um)
    F = 1; # (Hz)
    bl.mirror1.start_oscillation(axis='X', amp = A, frequency = F);
    bl.pause(15); # pauses after 10 seconds
    bl.mirror1.stop();
    
# Square Wave Generation
    # sampling_rate--> the number of point per second used to generate the square wave
    # sin_frequency --> the number of sin/square wave oscillations for a given period.
if False:
    sampling_rate = 100; # Hz --> it seems like a sample rate of 50 can capture enough points to construct a square wave.
    duration = 2; # seconds

    time = np.linspace(0, duration, int(duration * sampling_rate), endpoint=False);

    sin_frequency = 10; # frequency of the sin function.
    amp = 1; # amplitude of the square wave.
    waveform = amp*np.sign(np.sin(sin_frequency*time));

    plt.plot(time, waveform)
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title('Square Waveform')
    plt.show()
    
    function = waveform;
    num_points = len(waveform); # essentially equivalent to the frequency
 
# Triangle Wave Generation
    # sampling_rate--> the number of point per second used to generate the square wave
if False:
    sampling_rate = 50 # Hz
    duration = 2 # seconds

    time = np.linspace(0, duration, int(duration * sampling_rate), endpoint=False)

    triangle = 2*np.abs(2*time - 2*np.floor(time+0.5)) - 1

    plt.plot(time, triangle)
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title('Triangle Waveform')
    plt.show()
    
    function = triangle;
    num_points = len(triangle);
    
# Better Straigth Line Displacement Code
if False:
    start_position = mirror1.position
    amplitude = 2; # amplitude in microns
    period = 2; # define the period for periodic linear functions
    speed = amplitude/(period/4);
    number_periods = 2;
    #does the first amplitude displacement
    bl.mirror.move_by(dx = amplitude, dy = 0, speed = speed); # note that the first amplitude is set to be in the positive direction.
    
    # does all the displacements in between the first and last amplitude displacements
    while i < number_periods:
        bl.mirror1.move_by(dx = -2*amplitude, dy = 0, speed = speed);
        bl.mirror1.move_by(dx = 2*amplitude, dy = 0, speed = speed);
    
    # does the last amplitude displamcent
    bl.mirror.move_by(dx = -1*amplitude, dy = 0, speed = speed);
    end_position = mirror1.position
    print("The displacement is complete.");
    print("Trap 1s initial position is: ", start_position);
    print("Trap 1s final position is: ", end_position);
    
# The following loop takes the selected displacement function and determines the change in coordinates at each individual point used when plotting the function.
# len(D_x) should be equivalent to the sampling_rate * duration of the displacent.

D_x = []; # initiate the displacement vector
i = 0; #initiate the count variable

while i < (num_points - 1):

    D_x.append(function[i+1] - function[i]); 
    i = i + 1;

# - there may be some issue with passing long decimal values to the trap - the displacement vector has been rounded in the event that this poses an issues.
print(len(D_x)) # just to check that they have printed properly. 
#rounded_D_x = [round(val,6) for val in D_x]
#print(D_x)
#print(rounded_D_x)

# NOTE: The following code has been tested and it works - quantification of the ideal step_times and sampling freuquencies will be done by looking at the trap disaplements from the CTrap data.

# The following try and loop statements are used to pass the incremental displacements to the CTrap mirrors using internal bluelake functions.
# Pause times are included as to avoid all the displacements occuring extremely rapidly and "messing up" the desired waveform.
# The pause/step_time corresponds to the time delay between sampling points in the original function. 

step_time = 1/sampling_rate

# Trying to check the time it take to operate the for loop without the pause statement.

if False:
    
    try:
        start_time = time.monotonic();
        for dx in D_x:
            bl.pause(step_time); # Pause time before displacement for sin and triangle wave
            bl.mirror1.move.by(dx=dx, dy=0, speed = 0);     
            #bl.pause(step_time); # pause time after displacement for square wave.
    finally:
        bl.mirror1.stop();
        end_time = time.monotonic();      
        elapsed_time = end_time - start_time;
        estimated_operational_time = (elapsed_time)/num_points  # the total time it takes to operate the for loop without the pause statements divided by the number of operation in the loop = average time for a single operation. 
        print(estimated_operational_time);
        

# Ben G - 2023
# similar idea to pause time - however the bluklake pause function will not be used. Instead we will try to use the OS pause time.
# there may be a lower limit to the python.sleep function (on the order of 10 milisecond). This would allow for sampling on 100 Hz scale.

# NOTE: there seems to be a lower limit on the pause time. It might be hard to incrementally move the traps using a for loop at sampling frequencies that exceed 50 Hz. 
# the oscillation command works well for those higher order frequencies. 


if False: 
    
    try:
        
        for dx in D_X:
            time.sleep(0.01); # time is in seconds here
            bl.mirror1.move.by(dx=dx, dy=0, speed = 0);     
            
        finally:
            bl.mirror1.stop(); 
            print("Displacement complete");


# In[23]:


print('test');


# In[ ]:





# In[ ]:




