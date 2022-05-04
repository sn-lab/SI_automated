
%% imaging session setup
%set session name
baseDir = 'C:\Users\schafferlab\Desktop\Visual Cortex';
sessionName = 'test1';
objective_magnification = 25;

%connect to teensyvis
clear vs %clear "vis-struct"
% vs = teensyComm([], 'Find-Ports');
% assert(~isempty(vs.foundports),'Could not detect a connected teensyvis')
% vs.port = vs.foundports{1}
vs.port = 'COM7';
vs = teensyComm(vs, 'Connect');

%%
%(optional) align to a previously acquired template?
alignToPrevious = false;
previousTemplateName = '';

%create folder for new session
sessionDir = fullfile(baseDir,[datestr(now,'yyyy-mm-dd') ' ' sessionName]);
assert(~isfolder(sessionDir),['sessionName "' sessionDir '" already exists.'])
mkdir(sessionDir);

%%
%set typical imaging parameters
assert(strcmpi(hSI.acqState,'idle'),'scanimage is busy');   % make sure scanimage is in an idle state
hSI.hStackManager.enable = 0;
hSI.hRoiManager.scanZoomFactor = 1;     % define the zoom factor
hSI.hStackManager.framesPerSlice = 1;   % set number of frames to capture in one Grab
hSI.hStackManager.numSlices = 1;
hSI.hStackManager.numVolumes = 1;
hSI.hChannels.loggingEnable = 1;     % enable logging
hSI.hChannels.channelSave = [1 2 3 4];
hSI.hChannels.channelMergeColor = {'blue'  'green'  'gray'  'red'};
hSI.acqsPerLoop = 1;
wavelength = 0;
logFileCounter = 1;
channel_options.nch = 4; %number of channels in the source file
channel_options.chsh = [2 3 4]; %channels to use for registering shifts
channel_options.pr = 'max'; %projection type to use across channels ('max' or 'mean')
search_options.mpp = 43/objective_magnification; %microns per pixel (as scanimage understands it)
search_options.search_range = 50; %microns from the current location to search for the template
search_options.step_size = 5; %step size between slices
search_options.manual_check = 0; %whether to manually check the results


%% align imaging to previous template
if alignToPrevious
    [status, ~] = TemplateFinder(previousTemplateName,channel_options,search_options);
    assert(status,'template not found; exiting script');
end


%% take an image series at 1030 nm
wavelength = 0;
while wavelength~=1030
    answer = inputdlg('Enter current wavelength to continue:','Set laser source to satsuma 1030 nm',[1 80]);
    if ~isempty(answer)
        wavelength = str2double(answer{1});
    end
end

%set imaging parameters
color_duration = 10;
assert(strcmpi(hSI.acqState,'idle'),'scanimage is busy');   % make sure scanimage is in an idle state
hSI.hStackManager.enable = 0;
hSI.hScan2D.logFilePath = sessionDir;        % set the folder for logging Tiff files
hSI.hScan2D.logFileStem = '1030';      % set the base file name for the Tiff file
hSI.hScan2D.logFileCounter = logFileCounter;         % set the current Tiff file number
hSI.hRoiManager.linesPerFrame = 256;
hSI.hRoiManager.pixelsPerLine = 256;
fps = hSI.hRoiManager.scanFrameRate;
hSI.hStackManager.framesPerSlice = ceil(fps*color_duration);   % set number of frames to capture in one Grab
hSI.extTrigEnable = 0;

fprintf('Acquiring 1030 nm image series... ');
hSI.startGrab();                
while ~strcmpi(hSI.acqState,'idle') %could also be 'grab' (or 'loop'?)
    pause(1);
end
fprintf('done.\n');
pause(2);

% correct for motion and create a new template
imageName = fullfile(hSI.hScan2D.logFilePath,[hSI.hScan2D.logFileStem '_' sprintf(['%0' num2str(5) 'd'],logFileCounter) '.tif']);
channel_options.chsh = [2 3 4]; %channels to use for registering shifts
[status, templateName_1030, ~] = run_multichannel_normcorre('imageName',imageName,'channelOptions',channel_options);
assert(strcmp(status,'success'),'motion correction did not finish successfully');


%% take an image series at 920 nm
wavelength = 0;
while wavelength~=920
    answer = inputdlg('Enter current wavelength to continue:','Set laser source to chameleon 920 nm',[1 80]);
    if ~isempty(answer)
        wavelength = str2double(answer{1});
    end
end

% align to 1030 nm template
channel_options.chsh = [2 3];
search_options.search_range = 20;
search_options.step_size = 2;
[status, shifts] = TemplateFinder(templateName_1030,channel_options,search_options,hSI);
all_shifts(1,:) = shifts;
assert(status,'template not found; exiting script');
    
%set imaging parameters
color_duration = 10;
assert(strcmpi(hSI.acqState,'idle'),'scanimage is busy');   % make sure scanimage is in an idle state
hSI.hStackManager.enable = 0;
hSI.hScan2D.logFilePath = sessionDir;        % set the folder for logging Tiff files
hSI.hScan2D.logFileStem = '920';      % set the base file name for the Tiff file
hSI.hScan2D.logFileCounter = logFileCounter;         % set the current Tiff file number
hSI.hRoiManager.linesPerFrame = 256;
hSI.hRoiManager.pixelsPerLine = 256;
fps = hSI.hRoiManager.scanFrameRate;
hSI.hStackManager.framesPerSlice = ceil(fps*color_duration);   % set number of frames to capture in one Grab
hSI.extTrigEnable = 0;

fprintf('Acquiring 920 nm image series... ');
hSI.startGrab();                
while ~strcmpi(hSI.acqState,'idle') %could also be 'grab' (or 'loop'?)
    pause(1);
end
fprintf('done.\n');
pause(2);

% correct for motion and create a new template
imageName = fullfile(hSI.hScan2D.logFilePath,[hSI.hScan2D.logFileStem '_' sprintf(['%0' num2str(5) 'd'],logFileCounter) '.tif']);
channel_options.chsh = [2 3 4]; %channels to use for registering shifts
[status, templateName_920, ~] = run_multichannel_normcorre('imageName',imageName,'channelOptions',channel_options);
assert(strcmp(status,'success'),'motion correction did not finish successfully');


%% take a stack at 780
while wavelength~=780
    answer = inputdlg('Enter current wavelength to continue:','Set laser source to chameleon 780 nm',[1 80]);
    if ~isempty(answer)
        wavelength = str2double(answer{1});
    end
end

range = 100;
step_size = 1;
hSI.hScan2D.logFileStem = 'plaques';
hSI.hScan2D.logFileCounter = 1;
hSI.hStackManager.enable = 1;
hSI.hStackManager.centeredStack = 1;
hSI.hStackManager.stackDefinition = 'bounded';
hSI.hStackManager.stackMode = 'slow';
hSI.hStackManager.stackReturnHome = 1;
hSI.hStackManager.framesPerSlice = 1;
hSI.hStackManager.numVolumes = 1;
hSI.hRoiManager.linesPerFrame = 512;
hSI.hRoiManager.pixelsPerLine = 512;
hSI.hStackManager.boundedStackDefinition = 'stepSize';
hSI.hStackManager.stackZStepSize = step_size;
samplePosition = hSI.hMotors.samplePosition;

hSI.hMotors.moveSample(samplePosition + [0 0 -range]);
hSI.startFocus();
opts.Interpreter = 'tex';
opts.Default = 'Continue';
answer = input('Adjust laser power, then click "Set START". Press any key to continue');
if strcmpi(hSI.acqState,'idle')
    hSI.startFocus();
end
hSI.hMotors.moveSample(samplePosition + [0 0 2*range]);
answer = input('Adjust laser power, then click "Set STOP". Press any key to continue');
if ~strcmpi(hSI.acqState,'idle')
    answer = input('Abort current focus to continue. Press any key when done.');
end

fprintf('Acquiring 780 nm stack... ');
hSI.startGrab();                
while ~strcmpi(hSI.acqState,'idle') %could also be 'grab' (or 'loop'?)
    pause(1);
end
fprintf('done.\n');
pause(2);


%% take a stack at 920
while wavelength~=920
    answer = inputdlg('Enter current wavelength to continue:','Set laser source to chameleon 920 nm',[1 80]);
    if ~isempty(answer)
        wavelength = str2double(answer{1});
    end
end

range = 100;
step_size = 1;
hSI.hScan2D.logFileStem = 'vessels';
hSI.hScan2D.logFileCounter = 1;
hSI.hStackManager.enable = 1;
hSI.hStackManager.centeredStack = 1;
hSI.hStackManager.stackDefinition = 'bounded';
hSI.hStackManager.stackMode = 'slow';
hSI.hStackManager.stackReturnHome = 1;
hSI.hStackManager.framesPerSlice = 1;
hSI.hStackManager.numVolumes = 1;
hSI.hRoiManager.linesPerFrame = 512;
hSI.hRoiManager.pixelsPerLine = 512;
hSI.hStackManager.boundedStackDefinition = 'stepSize';
hSI.hStackManager.stackZStepSize = step_size;
samplePosition = hSI.hMotors.samplePosition;

hSI.hMotors.moveSample(samplePosition + [0 0 -range]);
hSI.startFocus();
opts.Interpreter = 'tex';
opts.Default = 'Continue';
answer = input('Adjust laser power, then click "Set START". Press any key to continue');
if strcmpi(hSI.acqState,'idle')
    hSI.startFocus();
end
hSI.hMotors.moveSample(samplePosition + [0 0 2*range]);
answer = input('Adjust laser power, then click "Set STOP". Press any key to continue');
if ~strcmpi(hSI.acqState,'idle')
    answer = input('Abort current focus to continue. Press any key when done.');
end

fprintf('Acquiring 920 nm stack... ');
hSI.startGrab();                
while ~strcmpi(hSI.acqState,'idle') %could also be 'grab' (or 'loop'?)
    pause(1);
end
fprintf('done.\n');
pause(2);


%% take a long image series at 920 (for spontaneous activity)
% align to 920 nm template
channel_options.chsh = [2 3 4];
search_options.search_range = 20;
search_options.step_size = 2;
[status, shifts] = TemplateFinder(templateName_920,channel_options,search_options);
all_shifts(2,:) = shifts;
assert(status,'template not found; exiting script');
    
%set imaging parameters
spontaneous_duration = 300; % seconds
assert(strcmpi(hSI.acqState,'idle'),'scanimage is busy');   % make sure scanimage is in an idle state
hSI.hStackManager.enable = 0;
hSI.hScan2D.logFilePath = sessionDir;        % set the folder for logging Tiff files
hSI.hScan2D.logFileStem = 'spont';      % set the base file name for the Tiff file
hSI.hScan2D.logFileCounter = logFileCounter;         % set the current Tiff file number
hSI.hRoiManager.linesPerFrame = 256;
hSI.hRoiManager.pixelsPerLine = 256;
fps = hSI.hRoiManager.scanFrameRate;
hSI.hStackManager.framesPerSlice = ceil(fps*spontaneous_duration);   % set number of frames to capture in one Grab
hSI.extTrigEnable = 0;

fprintf('Acquiring spontaneous image series... ');
hSI.startGrab();                
while ~strcmpi(hSI.acqState,'idle') %could also be 'grab' (or 'loop'?)
    pause(1);
end
fprintf('done.\n');
pause(2);

% correct for motion and create a new template
imageName = fullfile(hSI.hScan2D.logFilePath,[hSI.hScan2D.logFileStem '_' sprintf(['%0' num2str(5) 'd'],logFileCounter) '.tif']);
channel_options.chsh = [2 3 4]; %channels to use for registering shifts
[status, templateName_spont, ~] = run_multichannel_normcorre('imageName',imageName,'channelOptions',channel_options);
assert(strcmp(status,'success'),'motion correction did not finish successfully');


%% take a long, externally-triggered image series at 920 (for directional tuning)
% align to 920 nm template
% channel_options.chsh = [2 3 4];
% search_options.search_range = 20;
% search_options.step_size = 2;
% [status, shifts] = TemplateFinder(templateName_920,channel_options,search_options);
% all_shifts(3,:) = shifts;
% assert(status,'template not found; exiting script');

%set stimulus experiment parameters
vs.expname = 'stim';
vs.directory = sessionDir;
vs.trial_duration = 8;
vs.randomize = 1; %1=randomize order of conditions, 0=don't randomize
vs.num_trials = 13;
vs.num_reps = 1;
num_acqs = vs.num_trials*vs.num_reps;
angles = [-1 0:30:330]; %only angles are used in this script

% set starting/default grating parameters
default.patterntype = 1; %1 = square-gratings, 2 = sine-gratings, 3 = flicker
default.bar1color = [0 0 30]; %RGB color values of bar 1 [R=0-31, G=0-63, B=0-31]
default.bar2color = [0 0 0]; %RGB color values of bar 2
default.backgroundcolor = [0 0 15]; %RGB color values of background
default.barwidth = 20; % width of each bar (pixels) (1 pixel ~= 0.58 degrees)
default.numgratings = 4; % number of bright/dark bars in grating
default.angle = 0; % angle of grating (degrees) [0=drifting right, positive angles rotate clockwise]
default.frequency = 1.5; % temporal frequency of grating (Hz) [0.1-25]
default.position = [0, 0]; % x,y position of grating relative to display center (pixels)
default.predelay = 2; % delay after start command sent before grating pattern begins (s) [0.1-25.5]
default.duration = 2; % duration that grating pattern is shown (s) [0.1-25.5]
default.trigger = 0; % tells the teensy whether to wait for an input trigger signal (TTL) to start or not

%set imaging parameters
assert(strcmpi(hSI.acqState,'idle'),'scanimage is busy');   % make sure scanimage is in an idle state
hSI.hStackManager.enable = 0;
hSI.hScan2D.logFilePath = fullfile(sessionDir,'stim');        % set the folder for logging Tiff files
hSI.hScan2D.logFileStem = '';      % set the base file name for the Tiff file
hSI.hScan2D.logFileCounter = logFileCounter;         % set the current Tiff file number
hSI.hRoiManager.linesPerFrame = 256;
hSI.hRoiManager.pixelsPerLine = 256;
fps = hSI.hRoiManager.scanFrameRate;
hSI.hStackManager.framesPerSlice = ceil(fps*6);   % set number of frames to capture in one Grab
hSI.acqsPerLoop = num_acqs;
hSI.extTrigEnable = 1;
hSI.hScan2D.trigAcqInTerm = 'PFI0';
mkdir(hSI.hScan2D.logFilePath);

%create order in which conditions will be presented
if vs.randomize==1
    for r = 1:vs.num_reps
        order(r,:)=randperm(vs.num_trials);
    end
else
    order = repmat(1:vs.num_trials,[vs.num_reps 1]);
end
vs.order = order;

%start imaging
disp(['settings: ' num2str(vs.num_trials*vs.num_reps) ' trials, ' num2str(vs.trial_duration) ' s trial duration']);
fprintf('Acquiring direction tuning image series... ');
hSI.startLoop();    
pause(2)

%send stimuli
param = default;
for r = 1:vs.num_reps
    vs.rep = r; %keep this; helps data analaysis later
    for t = 1:vs.num_trials
        vs.trial = t; %keep this; helps data analaysis later
        if angles(order(r,t))==-1
            default.patterntype = 3;
            param.angle = 30*(r-1);
        else
            default.patterntype = 1;
            param.angle = angles(order(r,t));
        end
        
        tic
        vs = teensyComm(vs, 'Start-Pattern', param); %send pattern parameters and display the pattern
        while toc<vs.trial_duration %delay until next trial
            pause(0.001);
        end
        vs = teensyComm(vs, 'Get-Data'); %retrieve data sent from teensy about displayed pattern

    end
end
            
while ~strcmpi(hSI.acqState,'idle') %could also be 'grab' (or 'loop'?)
    pause(1);
end
fprintf('done.\n');
pause(2);

%save data for current experiment
filename = [datestr(now,'yyyy-mm-dd HH-MM-SS') ' ' vs.expname ' vs.mat'];
if ~exist(vs.directory,'dir')
    [pathstr,newfolder,~] = fileparts(vs.directory);
    mkdir(pathstr,newfolder);
end
save(fullfile(vs.directory,filename),'vs')

%close connection to display
vs = teensyComm(vs, 'Disconnect'); %close connection to controller

%%
%combine folder into a single file
[imageName, ~, ~] = combine_tifs({hSI.hScan2D.logFilePath}, [1 2 3 4]);

% correct for motion and create a new template
channel_options.chsh = [2 3 4]; %channels to use for registering shifts
[status, templateName_spont, ~] = run_multichannel_normcorre('imageName',imageName,'channelOptions',channel_options);
assert(strcmp(status,'success'),'motion correction did not finish successfully');

