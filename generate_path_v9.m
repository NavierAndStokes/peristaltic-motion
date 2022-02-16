%% Notes:
% This code is designed to calculate the path for the peristalsis
% experiment. The code produces the steps for a discretized motion in a text file.
% all units should be in cm/g/s
%% 
clear
close all

%% Details
sm = 1/4;                   % step mode 1- fullstep 1/2 - half step, and so on
thread = 1;                 % threads per inch of actuator
arduino_f = 1000;           % Arduino output frequency (takes 2 outputs to generate 1 step 010101)
ins = thread*200;           % threads per inch times number of steps
fr = 20*arduino_f/2;        % Number of signals that will be output seconds * the frequency of the output
maxs = arduino_f/ins/2;     % max speed of displacement 
dy = sm*1/(200*thread);     % step size

%% Desired path
% traveling wave example units are in cm
t = linspace(0,fr/(arduino_f/2),fr);
x = linspace(0,39.1160);
xd = linspace(0,39.1160,24);
yo = zeros(size(xd));

% Parameters
c = 25;                              % wave speed
d = 10;                              % some damping coefficient to initialize the motion
Amp = 2.54*0.2;                      % amplitude 
reseter = .95;                       % This allows some time for the motors to go back to the initial position to be reused .95 for 20 seconds seems to be good so 1 sec to recover approx
nu = 0.0089;                         % Viscosity
L = 2.54;                            % Channel gap
%% special

lambda = 30;
w = 2*pi*c/lambda;
lambda=c/(w/(2*pi));
alpha = 2*pi/lambda;
displacement_maxspeed = Amp*2*pi*c/lambda;
%% Estimated Reynolds

disp('estimated Reynolds number')
Re = L^2*w/nu    % viscosity of water  cm^2/s
disp('Wave parameters')
disp('Wave speed')
c
disp('Wavelength')
lambda
disp('w freq')
w
disp('Wave amplitude')
Amp
disp('Run time (s)')
max(t)

Parameters = [Re c lambda Amp max(t)];

%% Generate ideal path
disp('Generating wave...')
tic
% Preallocate memory
y = zeros(length(t),length(x));
yd = zeros(length(t),length(xd));

%% Sinewave
y = Amp*(1-exp(-t'/max(t)*d)).*(sin(2*pi/lambda*x+w*t'));  
yd = Amp*(1-exp(-t'/max(t)*d)).*(sin(2*pi/lambda*xd+w*t'));

% Stationary wave
% y = ones(size(x)).*Amp.*sin(w*(c*t'));
% yd = ones(size(xd)).*Amp.*sin(w*(c*t'));

% Move target only
% yd(:,21:24) = 0;
% yd(:,1:19) = 0;
%%
% reset to 0 point
y(round(reseter*length(t)):end,:) = 0;
yd(round(reseter*length(t)):end,:) = 0;



%% Determine correcting positon for discrete elements
% ideal position is y
% actual position is yo
% correcting position is yc
% Determine error in time step

yc = yo;
%% correct 
tolerance = dy;         % The tolerance generally is related to the step size
err_y = zeros(1,24);
Stepping = zeros(1,24);
yc = zeros(length(t),24);

for i=2:length(t)
    err_y = (yd(i-1,:)-yc(i-1,:));          % Difference between desired position (yd) and stepper position (yc)
    stepping = (abs(err_y) > tolerance);    % is the position too far from tolerance
    derr_y = err_y./abs(err_y);             % Direction of the step?
    derr_y(isnan(derr_y))=0;                % get rid of perfect positions?
    stepping = stepping.*derr_y;            % Combine logic and direction
    yc(i,:) = yc(i-1,:)+dy*stepping;        % New stepper position
end
%
toc
disp('Wave generation done')
disp('Plotting...')

%% Plot wave 
figure
for i=1:round(length(t)/2000):round(length(t))
    
   hold off
   plot(x,y(i,:),'linewidth',3) 
   hold on 
   ylim([-Amp Amp])
   plot(xd,yc(i,:),'rx','markersize',12,'linewidth',2) 
   plot(xd,yc(i,:),'k:','linewidth',2) 
   time = num2str(t(i),3);
   xlabel('cm')
   ylabel('cm')
   title(['t = ' time 's'])
   pause(.1)

end
disp('finished with generation of ideal wave')
%% convert into output
disp('generating outputable data ...')

dir = zeros(length(t),24);
step_m = zeros(length(t),24);

for i=2:length(t)
    err_y = (yd(i-1,:)-yc(i-1,:));
    for ei = 1:length(err_y)
        if abs(err_y(ei)) < tolerance
            err_y(ei) = 0;
        end
    end
    err_y = err_y./abs(err_y);
    err_y(isnan(err_y))=0;
    
    for ei = 1:length(err_y)
        if err_y(ei) ~= 0
            step_m(i,ei) = 1;
        end

        if err_y(ei) > 0
            dir(i,ei) = 1;
        end
    end
end

step_data = zeros(length(dir),48);
step_data(2:2:2*length(dir),25:25+23) = dir;
step_data(1:2:2*length(dir),25:25+23) = dir;
step_data(2:2:2*length(dir),1:24) = zeros(length(dir),24);
step_data(1:2:2*length(dir),1:24) = step_m;
disp('data done')

%% Saving outputable data and a parameter file
yon = input('save this? (1 = yes)');
if yon == 1
   
    namer = input('name?');
    csvwrite(['output_' namer '.txt'],step_data);
    csvwrite(['parameters_' namer '.txt'],Parameters);
    fid=fopen(['parameters_' namer '.txt'],'wt');
    fprintf(fid,'Re c lambda Amp max(t)');
    fprintf(fid,'\n');
    for i = 1:length(Parameters)
    fprintf(fid,' %f',Parameters(i));
    end
    fclose(fid);
end

%% Calculate volume change in the system due to wave

Yd = trapz(yd,2);
Yd = Yd*2.54; % area displaced*width of channel
figure(55)
plot(t,Yd)
title('Volume change in the system due to wave')
xlabel('time')
ylabel('cm^3')






