function [status, best_fit_shifts, search_options] = TemplateFinder(templateName,channel_options,search_options,hSI)
% FUNCTION TemplateFinder(templateName,channel_options,search_options)
%
% Automatically find and move to a location in the sample that most closely
% matches a supplied template image.
% 
% INPUTS
% templateName: full filename of the template .tif file
% channel_options.nch: number of channels in the source file
% channel_options.chsh: channels to use for registering shifts
% channel_options.pr: projection type to use across channels ('max' or 'mean')
% search_options.search_range: microns from the current location to search for the template
% search_options.step_size: step size between slices
% search_options.fit_method: 'poly-fit' or 'max'
% search_options.manual_check: whether to manually check the results
%
% OUTPUTS
% status: outcome of the templatefinder function
% best_fit_shifts: distance between starting location and template location (if found)


%set parameters
assert(strcmpi(hSI.acqState,'idle'),'scanimage is busy');   % make sure scanimage is in an idle state
zPos = -search_options.search_range:search_options.step_size:search_options.search_range;
num_slices = length(zPos);
samplePosition = hSI.hMotors.samplePosition; %get current sample position
framesPerSlice = 5;
numChannels = 4;
[optimizer, metric] = imregconfig('multimodal');

%load the template
template = read_file(templateName);
[templateHeight,templateWidth] = size(template);

%% if ppm isn't supplied, calculate it here:
if ~isfield(search_options,'ppm') || isempty(search_options.ppm)
    %set imaging parameters for single images
    assert(strcmpi(hSI.acqState,'idle'),'scanimage is busy');   % make sure scanimage is in an idle state
    hSI.hStackManager.enable = 0;
    hSI.hRoiManager.scanZoomFactor = 1;     % define the zoom factor
    hSI.hStackManager.framesPerSlice = 1;   % set number of frames to capture in one Grab
    hSI.hStackManager.numSlices = 1;
    hSI.hStackManager.numVolumes = 1;
    hSI.hChannels.loggingEnable = 1;     % enable logging
    hSI.hChannels.channelSave = [1 2 3 4];
    hSI.hScan2D.logFileStem = 'tempGrab';      % set the base file name for the Tiff file
    hSI.hScan2D.logFileCounter = 1;         % set the current Tiff file number
    hSI.hRoiManager.linesPerFrame = templateHeight;
    hSI.hRoiManager.pixelsPerLine = templateWidth;
    hSI.extTrigEnable = 0;
    moveDist = round(0.2*min([templateHeight templateWidth]));

    %acquire images at 3 nearby locations at the same depth
    fprintf('Acquiring images to calibrate pixels/micron... ');
    hSI.startGrab(); %starting location          
    while ~strcmpi(hSI.acqState,'idle') %could also be 'grab' (or 'loop'?)
        pause(1);
    end
    hSI.hScan2D.logFileCounter = 2;
    hSI.hMotors.moveSample(samplePosition + [moveDist 0 0]);
    hSI.startGrab(); %moved small amount in x direction                
    while ~strcmpi(hSI.acqState,'idle') %could also be 'grab' (or 'loop'?)
        pause(1);
    end
    hSI.hScan2D.logFileCounter = 3;
    hSI.hMotors.moveSample(samplePosition + [0 moveDist 0]);
    hSI.startGrab(); %moved small amount in y direction            
    while ~strcmpi(hSI.acqState,'idle') %could also be 'grab' (or 'loop'?)
        pause(1);
    end
    hSI.hMotors.moveSample(samplePosition);  %move back to start 
    fprintf('done.\n');

    %load the 3 images
    grabs = nan([templateHeight templateWidth 3]);
    for i = 1:3
        grabName = fullfile(hSI.hScan2D.logFilePath,[hSI.hScan2D.logFileStem '_' sprintf(['%0' num2str(5) 'd'],i) '.tif']);
        grab = read_file(grabName);
        grabs(:,:,i) = max(grab(:,:,channel_options.chsh),[],3,'omitnan'); %create composite frame
        delete(grabName)
    end

    %calculate the pixels-per-micron conversion factors for the current state of the microscope
    %a*x = b
    %x = a\b - matrix left division to solve system of linear equations
    %x = micron movements (e.g. [1;1])
    %a = ppm (e.g. [ppm_x' ppm_y'])
    %b = pixel movements (e.g. imr_tform.T(3,[1 2])')
    imr_tform = imregtform(grabs(:,:,2),grabs(:,:,1),'translation',optimizer,metric); %movement in pixels
    search_options.ppm(:,1) = imr_tform.T(3,[1 2])'/moveDist; %ppm (for x micron movements)
    imr_tform = imregtform(grabs(:,:,3),grabs(:,:,1),'translation',optimizer,metric);
    search_options.ppm(:,2) = imr_tform.T(3,[1 2])'/moveDist; %ppm (for y micron movements)
end


%set imaging parameters for a stack
hSI.hScan2D.logFileStem = 'tempstack';
hSI.hScan2D.logFileCounter = 1;
hSI.hStackManager.enable = 1;
hSI.hStackManager.centeredStack = 1;
hSI.hStackManager.framesPerSlice = framesPerSlice;
hSI.hStackManager.stackDefinition = 'uniform';
hSI.hStackManager.stackMode = 'slow';
hSI.hStackManager.stackReturnHome = 1;
hSI.hStackManager.stackZStepSize = search_options.step_size;
hSI.hStackManager.numSlices = num_slices;
shifts = nan(num_slices,2); %shifts in [x, y, z] directions
stackName = fullfile(hSI.hScan2D.logFilePath,[hSI.hScan2D.logFileStem '_' sprintf(['%0' num2str(5) 'd'],1) '.tif']);
if isfile(stackName)
    delete(stackName)
end

%%
%acquire a stack
hSI.startGrab();      
fprintf('Acquiring a stack centered on the current depth... ');
while ~strcmpi(hSI.acqState,'idle')
    pause(1);
end
fprintf('done.\n');

%load the stack 
tempstack = read_file(stackName);
tempstack = reshape(tempstack,[templateHeight templateWidth numChannels framesPerSlice num_slices]);
tempstack = mean(tempstack,4); %average framesPerSlice (if >1)
tempstack = max(tempstack(:,:,channel_options.chsh,:,:),[],3,'omitnan');  %create composite frames
tempstack = permute(tempstack,[1 2 5 3 4]);

%register each frame of the stack to the template
R = nan(1,num_slices);
str = num2str(1);
fprintf(['Registering ' num2str(num_slices) ' slices to the template... ' str])
for s = 1:num_slices
    fprintf([repmat('\b',[1 length(str)]) num2str(s)])
    str = num2str(s);
    frame = tempstack(:,:,s);
    
    %perform rigid motion registration on the composite to the template
    imr_tform = imregtform(frame,template,'translation',optimizer,metric);
    frame_reg = imwarp(frame,imr_tform,'OutputView',imref2d(size(template)));

    %get correlation between registered template2 and template1
    kept_mask = imwarp(ones(size(template)),imr_tform,'OutputView',imref2d(size(template)));
    template_masked = template; 
    template_masked(~kept_mask) = 0;
    R(s) = corr2(template_masked,frame_reg);
    
    %calculate x, y, and rotation shifts
    shifts(s,:) = imr_tform.T(3,[1 2]); 
end
fprintf([repmat('\b',[1 length(str)]) 'done.\n'])

%delete tempstack file
if isfile(stackName)
    delete(stackName)
end


%% find frame with max correlation coefficients
switch search_options.fit_method
%     case 'poly-fit'
%         %fit polynomial to correlation-depth data to find the best depth
%         z = shifts(:,3)';
%         fitobject = fit(z',R','poly2');
%         coeffvals= coeffvalues(fitobject);
%         highresZ = -search_range:search_range;
%         highresR = coeffvals(1)*(highresZ).^2 + coeffvals(2)*(highresZ) + coeffvals(3);
%         [max_R,~] = max(highresR);
%         
%         
%         hSI.hMotors.moveSample(samplePosition + [0 0 max_R]); %move to the best depth
%         
%         %%%%take an image
%         
%         hSI.hMotors.moveSample(samplePosition + [0 0 -max_R]); %move back
%         %register image to template
%         %calculate x,y shifts
%         %save best_frame_reg

    case 'max'
        [max_R,max_R_ind] = max(R);
        best_fit_shift_pixels = shifts(max_R_ind,:)';
        frame = tempstack(:,:,max_R_ind); 
        imr_tform = imregtform(frame,template,'translation',optimizer,metric);
        best_frame_reg = double(imwarp(frame,imr_tform,'OutputView',imref2d(size(template))));
        
    otherwise
        error('Only ''max'' fit method is implemented so far') 
end
%%
%convert shifts from pixels to x and y movements
best_fit_shifts = -round(search_options.ppm\best_fit_shift_pixels)';
best_fit_shifts(3) = zPos(max_R_ind);
shift_string = sprintf('shifts: \nx=%d, y=%d, z=%d\n',best_fit_shifts);
    
if max_R>0.95
    fprintf(['template location found! (R = ' num2str(max_R) ')\n']);
    fprintf(shift_string);
    status = 1;
elseif max_R>0.8
    fprintf(['possible template location found. (R = ' num2str(max_R) ')\n']);
    fprintf(shift_string);
    status = 2;
else
    fprintf(['low confidence that template location was found. (max R = ' num2str(max_R) ')\n']);
    status = 0;
end
    
if search_options.manual_check
    manualFig = figure();
    imshow([template/max(template,[],'all') best_frame_reg/max(best_frame_reg,[],'all')])
    %ask for input to decide status
    answer = input('Continue with current result, or manually adjust? (1=continue, 2=manual, 3=stop script): ');
    switch answer
        case 1
            status = 1;
        case 2
            answer = input('Try manually adjusting, then decide the next step. (1=continue, 3=stop script): ');
            switch answer
                case 1
                    status = 3;
                case 3
                    status = 0;
                otherwise
                    error('Only 1 and 3 are valid options')
            end
        case 3
            status = 0;
        otherwise
            error('Only 1, 2, and 3 are valid options')
    end
    if ishghandle(manualFig)
        close(manualFig);
    end
end

%% move the motors based on the shifts
if status
    hSI.hMotors.moveSample(samplePosition + best_fit_shifts); 
    fprintf('Repositioned sample to match the template location.\n');
end
