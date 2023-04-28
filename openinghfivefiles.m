% Ben G. + Derek A. 2023
% This code takes a LUMICKS h5 files and converts the requires inputs into
% a 4 column txt file (x_QUAD, y_QUAD, x_AOD, y_AOD). This format fits the current format being used.
% LUMICKS data from the h5 file comes out dirrectly as force and trap
% position. The trap response (pN/V) is used to convert the the force
% measurement back into a voltage that can then be used to create a power
% spectrum. 
close all
clear all
clc
numBeads = input('how many traps did you use to trap beads (1 - one trap, 2 - two traps)? ');

file = '4-14/20230415-172622 Marker 3.h5'; 
h5disp(file);

%extracting force, position, time and sampling data from the h5 files obtained via the
%ctrap bluelake software
force_1x = h5read(file, '/Force HF/Force 1x');
force_1y = h5read(file, '/Force HF/Force 1y');


trap_position_1x = h5read(file, '/Trap position/1X');
trap_position_1y = h5read(file, '/Trap position/1Y');
offset1x = h5readatt(file, '/Calibration/1/Force 1x','Offset (pN)');
response1x = h5readatt(file, '/Calibration/1/Force 1x','Response (pN/V)'); % this is going to be used to get the volts
offset1y = h5readatt(file, '/Calibration/1/Force 1y','Offset (pN)');
response1y = h5readatt(file, '/Calibration/1/Force 1y','Response (pN/V)');
volts_1x = (force_1x -offset1x)./response1x;
volts_1y = (force_1y - offset1y)./response1y;
x1AOD = trap_position_1x;% trap position (um)
y1AOD = trap_position_1y;% trap position (um)

if numBeads == 2
    force_2x = h5read(file,'/Force HF/Force 2x');
    force_2y = h5read(file,'/Force HF/Force 2y');
    %c = trap_position_1x(1); %track initial coordinates
    c = 30.22;
    d = trap_position_1y(5000000);
    trap_position_2x = repelem(c,length(trap_position_1x)); %this is dependent on the trial you are running
    trap_position_2x = trap_position_2x.';
    trap_position_2y = repelem(d-3,length(trap_position_1y)); %this is also dependent on trial you are running
    trap_position_2y = trap_position_2y.';
    offset2x = h5readatt(file, '/Calibration/2/Force 2x','Offset (pN)');
    response2x = h5readatt(file, '/Calibration/2/Force 2x','Response (pN/V)'); % this is going to be used to get the volts
    offset2y = h5readatt(file, '/Calibration/2/Force 2y','Offset (pN)');
    response2y = h5readatt(file, '/Calibration/2/Force 2y','Response (pN/V)');
    volts_2x = (force_2x -offset2x)./response2x;
    volts_2y = (force_2y - offset2y)./response2y;
    x2AOD = trap_position_2x;% trap position (um)
    y2AOD = trap_position_2y;% trap position (um)
end




%The lumicks workflow requires us to determine the start and stop times
%when creating a power spectrum. Adam's workflow does not require these
%inputs and are thus commented out. The following code gets the start and
%stop time of the run (in ns), the number of individual points, and the
%timestep of each sampling point. 

%start = h5readatt(file,'/Force HF/Force 1x','Start time (ns)');
%stop = h5readatt(file,'/Force HF/Force 1x','Stop time (ns)');
%num_points = length(force_data_1x);
%sampling_rate = h5readatt(file,'/Force HF/Force 1x','Sample rate (Hz)');
%time_step = 1/sampling_rate;
%t = [0:time_step:0+(num_points-1)*time_step];

% initializes the variable arrays that will represent the columns in the
% final txt file
%volts_x = (force_1x-offsetx)/responsex; %xQUAD (V)
%volts_y = (force_1y-offsety)/responsey; %yQUAD (V)

if numBeads == 2
    
    raw_data_table = [volts_1x,volts_1y, x1AOD, y1AOD,volts_2x,volts_2y,x2AOD,y2AOD]; % variable are entered into a matrix prior to be written into a txt file.
    writematrix(raw_data_table,'3 um sep y.txt','Delimiter','tab') 
elseif numBeads == 1
    raw_data_table = [volts_1x,volts_1y, x1AOD, y1AOD]; % variable are entered into a matrix prior to be written into a txt file.
    writematrix(raw_data_table,'multifrequency bead oscillation 428.txt','Delimiter','tab') 
end




%figure
hold off
positionBead = ((volts_1x *1496)./1000);
positionBeadrelative = positionBead- mean(positionBead);
positionBead = positionBeadrelative + trap_position_1x;
t = 1:8546575;
figure
plot(t,positionBeadrelative)
axis([0 5000000 -0.05 0.05])
xlabel('sample index @ 78125 Hz')
ylabel('position of bead from centerpoint (um)')
hold off
figure
plot(t(1:5000000),trap_position_1x(1:5000000))
xlabel('sample index @ 78125 Hz')
ylabel('global position of trap (um)')


%{
plot(a, trap_position_1x + position)
ylabel("bead position (microns)")
%}
%}