

%% take a bounded stack
hSI.hScan2D.logFileStem = 'stack';
hSI.hStackManager.enable = 1;
hSI.hStackManager.centeredStack = 0;
hSI.hStackManager.stackDefinition = 'bounded';
hSI.hStackManager.stackMode = 'slow';
hSI.hStackManager.useStartEndPowers = 1;
hSI.hStackManager.stackReturnHome = 1;
% hSI.hStackManager.stackStartPowerFraction = [0.05 0]; %read only
% hSI.hStackManager.stackEndPowerFraction = [0.1 0];
hSI.hStackManager.framesPerSlice = 1;
hSI.hStackManager.numSlices = 3;
hSI.hStackManager.numVolumes = 1;
% hSI.hStackManager.stackZEndPos = 10; %have to be set in the gui?
% hSI.hStackManager.stackZStartPos = -10;

fprintf('Acquiring... ');

hSI.startGrab(); 
while ~strcmpi(hSI.acqState,'idle')
    pause(1);
end
fprintf('done.\n');


%% Take multiple externally-triggered image series'
hSI.hStackManager.enable = 0;
% hSI.hMotors.samplePosition = [0 0 0];    % move stage to origin Note: depending on motor this value is a 1x3 OR 1x4 matrix
hSI.hScan2D.logFilePath = 'C:\Users\schafferlab\Desktop\Visual Cortex';        % set the folder for logging Tiff files
hSI.hScan2D.logFileStem = 't-iseries';      % set the base file name for the Tiff file
hSI.hScan2D.logFileCounter = 1;         % set the current Tiff file number
hSI.hRoiManager.scanZoomFactor = 1.5;     % define the zoom factor
hSI.hStackManager.framesPerSlice = 10;   % set number of frames to capture in one Grab
hSI.hStackManager.numSlices = 1;
hSI.hStackManager.numVolumes = 1;
hSI.acqsPerLoop = 3;
hSI.hRoiManager.linesPerFrame = 256;
hSI.hRoiManager.pixelsPerLine = 256;
% hSI.hRoiManager.scanFrameRate
hSI.hChannels.loggingEnable = 1;     % enable logging
hSI.hChannels.channelSave = [1 2 3 4];
hSI.hChannels.channelMergeColor = {'blue'  'green'  'gray'  'red'};
hSI.extTrigEnable = 1;
hSI.hScan2D.trigAcqInTerm = 'PFI0';

fprintf('Acquiring... ');

hSI.startGrab();                        
% hSI.startLoop(); %maybe?
while ~strcmpi(hSI.acqState,'idle') %could also be 'grab' (or 'loop'?)
    pause(1);
end
fprintf('done.\n');


%% take an image series, correct for motion, create a new template
assert(strcmpi(hSI.acqState,'idle'),'scanimage is busy');   % make sure scanimage is in an idle state
hSI.hStackManager.enable = 0;
% hSI.hMotors.samplePosition = [0 0 0];    % move stage to origin Note: depending on motor this value is a 1x3 OR 1x4 matrix
hSI.hScan2D.logFilePath = 'C:\Users\schafferlab\Desktop\Visual Cortex';        % set the folder for logging Tiff files
hSI.hScan2D.logFileStem = 'iseries';      % set the base file name for the Tiff file
hSI.hScan2D.logFileCounter = 1;         % set the current Tiff file number
hSI.hRoiManager.scanZoomFactor = 1.5;     % define the zoom factor
hSI.hStackManager.framesPerSlice = 20;   % set number of frames to capture in one Grab
hSI.hStackManager.numSlices = 1;
hSI.hStackManager.numVolumes = 1;
hSI.acqsPerLoop = 1;
hSI.hRoiManager.linesPerFrame = 256;
hSI.hRoiManager.pixelsPerLine = 256;
% hSI.hRoiManager.scanFrameRate
hSI.hChannels.loggingEnable = 1;     % enable logging
hSI.hChannels.channelSave = [1 2 3 4];
hSI.hChannels.channelMergeColor = {'blue'  'green'  'gray'  'red'};
hSI.extTrigEnable = 0;

fprintf('Acquiring... ');
hSI.startGrab();                        
% hSI.startLoop(); %maybe?
while ~strcmpi(hSI.acqState,'idle') %could also be 'grab' (or 'loop'?)
    pause(1);
end
fprintf('done.\n');

imageName = fullfile(hSI.hScan2D.logFilePath,[sprintf('%04d',hSI.hScan2D.logFileStem) '.tif']); %%%%%%%%%%%%%%check this
%settings for current image
channel_options.nch = 4; %number of channels in the source file
channel_options.chsh = [2 3 4]; %channels to use for registering shifts
channel_options.pr = 'max'; %projection type to use across channels ('max' or 'mean')

%perform multichannel normcorre motion correction and create a new template
[status, templateName, ~] = run_multichannel_normcorre('imageName',imageName,...
    'channelOptions',channel_options);


%% take a stack; register each frame to a template; find the best location
assert(strcmpi(hSI.acqState,'idle'),'scanimage is busy');   % make sure scanimage is in an idle state
search_range = 50; %microns from the current location to search for the template
step_size = 5;
num_channels = 4;
comp_channels = [2 3 4];
samplePosition = hSI.hMotors.samplePosition; %get current position
axesPosition = hSI.hMotors.axesPosition;

%load the template
template = read_file(templateName);
[templateHeight,templateWidth] = size(template);

hSI.hScan2D.logFileStem = 'tempstack';
hSI.hScan2D.logFileCounter = 1;
hSI.hStackManager.enable = 1;
hSI.hStackManager.centeredStack = 1;
hSI.hStackManager.stackDefinition = '___________'; %the 1st option
hSI.hStackManager.stackMode = 'slow';
hSI.hStackManager.stackReturnHome = 1;
hSI.hStackManager.framesPerSlice = 1;
hSI.hRoiManager.linesPerFrame = templateHeight;
hSI.hRoiManager.pixelsPerLine = templateWidth;
%%%%%%%%set start, stop, step size
hSI.hStackManager.numVolumes = 1;
numSlices = hSI.hStackManager.numSlices; %%%%%%%%%%
shifts = nan(numSlices,4); %shifts in [x, y, z, rotation] directions
shifts(:,3) = -search_range:step_size:search_range; %%%%%%%%%%%%%%%%% z-pos from start to end

hSI.startGrab();                        
% hSI.startLoop(); %maybe?
while ~strcmpi(hSI.acqState,'idle') %could also be 'grab' (or 'loop'?)
    pause(1);
end
fprintf('done.\n');

%load tempstack
stackName = fullfile(hSI.hScan2D.logFilePath,[sprintf('%04d',hSI.hScan2D.logFileStem) '.tif']); %%%%%%%%%%%%%%check this
tempstack = read_file(fullfile(stackName));

%register each frame of the stack to the template
[optimizer, metric] = imregconfig('multimodal');
R = nan(1,numSlices);
str = num2str(1);
fprintf(['Registering ' num2str(numSlices) ' frames to the template... ' str])
for s = 1:numSlices
    fprintf([repmat('\b',[1 length(str)]) num2str(s)])
    str = num2str(s);
    frames = comp_channels + num_channels*(s-1);
    frame = max(tempstack(:,:,frames),[],3,'omitnan'); %create composite frame
    
    %perform rigid motion registration on the composite to the template
    imr_tform = imregtform(frame,template,'translation',optimizer,metric);
    frame_reg = imwarp(frame,imr_tform,'OutputView',imref2d(size(template)));

    %get correlation between registered template2 and template1
    kept_mask = imwarp(ones(size(template)),imr_tform,'OutputView',imref2d(size(template)));
    template_masked = template; 
    template_masked(~kept_mask) = 0;
    R(s) = corr2(template_masked,frame_reg);
    
    %calculate x, y, and rotation shifts
    shifts(s,1) = imr_tform.T(3,1); 
    shifts(s,2) = imr_tform.T(3,2);
    shifts(s,4) = 0; %%%%%%%%%%rotation in degrees (switch to rigid later?)
end
fprintf([repmat('\b',[1 length(str)]) 'done.\n'])

%find frame with max correlation coefficients
z = shifts(:,3)';
highresZ = -search_range:step_size/10:search_range;
highresR = coeffvals(1)*(highresZ).^2 + coeffvals(2)*(highresZ) + coeffvals(3);
fitobject = fit(z',R','poly2');
coeffvals= coeffvalues(fitobject);
fit_R = coeffvals(1)*(z).^2 + coeffvals(2)*(z) + coeffvals(3);
[~,max_R_ind] = max(fit_R);
max_R = R(max_R_ind); 
best_fit_shifts = shifts(max_R_ind,:);
shift_string = sprintf('shifts:\nx=%d, y=%d, z=%d, rot=%d\n',best_fit_shifts);

if max_R>0.95
    fprintf(['template location found! (R=' num2str(max_R) ')\n']);
    fprintf(shift_string);
    if abs(best_fit_shifts(4))>10
        warning(['location appears to be rotated relative to template; rotate the sample ' num2str(-best_fit_shifts(4)) ' deg CW to fix.'])
    end
elseif max_R>0.85
    fprintf(['possible template location found. (R=' num2str(max_R) ')\n']);
    fprintf(shift_string);
    if abs(best_fit_shifts(4))>10
        warning(['location appears to be rotated relative to template; rotate the sample ' num2str(-best_fit_shifts(4)) ' deg CW to fix.'])
    end
else
    fprintf(['template location not found. (max R=' num2str(max_R) ')\n']);
    if abs(best_fit_shifts(4))>10
        warning(['location might be rotated relative to template; rotate the sample ' num2str(-best_fit_shifts(4)) ' deg CW to fix.'])
    end
end


%% move the motors based on the shifts (write message for rotation?)
hSI.hMotors.samplePosition = samplePosition + best_fit_shifts(1:3);
fprintf('Repositioned sample to match the template location.\n');

