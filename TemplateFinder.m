function [status, best_fit_shifts] = TemplateFinder(templateName,channel_options,search_options,hSI)
% FUNCTION TemplateFinder(templateName,channel_options,search_options)
% take a stack; register each frame to a template; find the best location
%%%%%%%%%%%%%make sure to delete the old temp stuff
% search_options.mpp = 2.15;
% search_options.search_range = 50; %microns from the current location to search for the template
% search_options.step_size = 5; %step size between slices
% search_options.manual_check = 0; %whether to manually check the results
% channel_options.nch = 4; %number of channels in the source file
% channel_options.chsh = [2 3 4]; %channels to use for registering shifts
% channel_options.pr = 'max'; %projection type to use across channels ('max' or 'mean')

%set parameters
assert(strcmpi(hSI.acqState,'idle'),'scanimage is busy');   % make sure scanimage is in an idle state
zPos = -search_options.search_range:search_options.step_size:search_options.search_range;
num_slices = length(zPos);
samplePosition = hSI.hMotors.samplePosition; %get current sample position

%load the template
template = read_file(templateName);
[templateHeight,templateWidth] = size(template);

hSI.hScan2D.logFileStem = 'tempstack';
hSI.hScan2D.logFileCounter = 1;
hSI.hStackManager.enable = 1;
hSI.hStackManager.centeredStack = 1;
hSI.hStackManager.stackDefinition = 'uniform';
hSI.hStackManager.stackMode = 'slow';
hSI.hStackManager.stackReturnHome = 1;
hSI.hStackManager.framesPerSlice = 1;
hSI.hStackManager.numVolumes = 1;
hSI.hRoiManager.linesPerFrame = templateHeight;
hSI.hRoiManager.pixelsPerLine = templateWidth;
hSI.hStackManager.stackZStepSize = search_options.step_size;
hSI.hStackManager.numSlices = num_slices;
shifts = nan(num_slices,4); %shifts in [x, y, z, rotation] directions
shifts(:,3) = zPos; % z-pos from start to end
stackName = fullfile(hSI.hScan2D.logFilePath,[hSI.hScan2D.logFileStem '_' sprintf(['%0' num2str(5) 'd'],1) '.tif']);
if isfile(stackName)
    delete(stackName)
end

hSI.startGrab();      
fprintf('Acquiring... ');
while ~strcmpi(hSI.acqState,'idle')
    pause(1);
end
fprintf('done.\n');

%load tempstack 
tempstack = read_file(fullfile(stackName));

%register each frame of the stack to the template
[optimizer, metric] = imregconfig('multimodal');
R = nan(1,num_slices);
str = num2str(1);
fprintf(['Registering ' num2str(num_slices) ' frames to the template... ' str])
for s = 1:num_slices
    fprintf([repmat('\b',[1 length(str)]) num2str(s)])
    str = num2str(s);
    frames = channel_options.chsh + channel_options.nch*(s-1);
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
    shifts(s,4) = 0; %rotation in degrees (switch to rigid registration later?)
end
fprintf([repmat('\b',[1 length(str)]) 'done.\n'])

%delete tempstack file
if isfile(stackName)
    delete(stackName)
end

%% find frame with max correlation coefficients
% fit method
% z = shifts(:,3)';
% fitobject = fit(z',R','poly2');
% coeffvals= coeffvalues(fitobject);
% highresZ = -search_range:step_size/10:search_range;
% highresR = coeffvals(1)*(highresZ).^2 + coeffvals(2)*(highresZ) + coeffvals(3);
% fit_R = coeffvals(1)*(z).^2 + coeffvals(2)*(z) + coeffvals(3);
% [~,max_R_ind] = max(fit_R);
% max_R = R(max_R_ind); 

% max method
[max_R,max_R_ind] = max(R);
best_fit_shifts = shifts(max_R_ind,:);
shift_string = sprintf('shifts:\nx=%d, y=%d, z=%d, rot=%d\n',best_fit_shifts);
    
if search_options.manual_check
    s = max_R_ind;
    frames = channel_options.chsh + channel_options.nch*(s-1);
    frame = max(tempstack(:,:,frames),[],3,'omitnan'); %create composite frame
    imr_tform = imregtform(frame,template,'translation',optimizer,metric);
    frame_reg = imwarp(frame,imr_tform,'OutputView',imref2d(size(template)));
    imshow([template/max(template,'all') frame_reg/max(frame_reg,'all')])
    %ask for input to decide status
    opts.Interpreter = 'tex';
    opts.Default = 'Continue';
    answer = questdlg('Do you want to continue with the TemplateFinder result, or would you like to find the template manually?', ...
	'Options to proceed', ...
	'Continue','Manual adjustment','Stop script',opts);
    switch answer
        case 'Continue'
            status = 1;
        case 'Manual adjustment'
            answer = questdlg('Manually adjust, then continue when finished', ...
                'Waiting for manual adjustment', ...
                'Continue','Stop script',opts);
            switch answer
                case 'Continue'
                    status = 3;
                case 'Stop script'
                    status = 0;
            end
        case 'Stop script'
            status = 0;
    end
else
    if max_R>0.95
        fprintf(['template location found! (R=' num2str(max_R) ')\n']);
        fprintf(shift_string);
        status = 1;
        if abs(best_fit_shifts(4))>10
            warning(['location appears to be rotated relative to template; rotate the sample ' num2str(-best_fit_shifts(4)) ' deg CW to fix.'])
        end
    elseif max_R>0.8
        fprintf(['possible template location found. (R=' num2str(max_R) ')\n']);
        fprintf(shift_string);
        status = 2;
        if abs(best_fit_shifts(4))>10
            warning(['location appears to be rotated relative to template; rotate the sample ' num2str(-best_fit_shifts(4)) ' deg CW to fix.'])
        end
    else
        fprintf(['template location not found. (max R=' num2str(max_R) ')\n']);
        status = 0;
        if abs(best_fit_shifts(4))>10
            warning(['location might be rotated relative to template; rotate the sample ' num2str(-best_fit_shifts(4)) ' deg CW to fix.'])
        end
    end
end

%% move the motors based on the shifts (write message for rotation?)
if status
    hSI.hMotors.moveSample(samplePosition + [search_options.mpp*best_fit_shifts([1 2]) best_fit_shifts(3)]); 
    fprintf('Repositioned sample to match the template location.\n');
end
