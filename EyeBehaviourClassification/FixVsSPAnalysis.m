function [classification] = FixVsSPAnalysis(SacStructInput)
% This script takes the processed Eye movemnt data from SaccadeProcessing, and
% classifies the foveations as either fixations or smooth pursuits. This is
% modified from the original human eye movement processing presented in
% Larsson et al. (2015).

% The script 
% 1) Remove any fast data not caught by the saccade processing script
% 2) find the beginning and end of each foveation, remove ones too short to
%    process
% 3) Measure the directionality of the foveations, and separate into
%    sections based on presence or absence of directionality (or whole 
%    section if no consistent directionality)
% 4) Calculate Dispersion, Consistency, Path ratio and Spatial Range for section
% 5) Classify section as Fixation, SP or leave undetermined
% 6) Process the undetermined sections by grouping adjacent undetermined
%    sections 
% 7) Calculate mean direction of all sections, and group with similar sections
% 8) Classify based on Path Ratio and Spatial Range


%% initialize the variable to hold the final classification
classification = zeros(length(SacStructInput),1);

%thresholds for various tests
DispersionThreshold = .45;
ConsistDirThreshold = .50;
PathRatioThreshold = .30;
SpatialRangeThreshold = 1.5; % max spatial range in degrees for fixation
MinSpatialRange = 1.0; % in degrees for smooth pursuit

%% 1) Remove any fast data not caught by the saccade processing script

velocity = medfilt1(sqrt((diff(SacStructInput(:,1)).^2)+(diff(SacStructInput(:,2)).^2))./...
    diff(SacStructInput(:,3)));
Fast = velocity>100;
% NaN everything above 100 degrees/s, something got missed
SacStructInput(Fast,1:2) = NaN;

%% 2) find the beginning and end of each foveation, 
%break up all intersaccadic intervals and analyse individually
goodData = [0;~isnan(SacStructInput(:,1));0];
startInds = find(diff(goodData)==1); %finds the start of every section of analyzable data
endInds = find(diff(goodData)==-1)-1; %finds the end of every section of analyzable data
%remove sections taht are too short
AllInds = [startInds endInds];
AllInds(abs(diff(AllInds,1,2))<20,:)=[];

startInds = AllInds(:,1);
endInds = AllInds(:,2);

%% 3) Process each foveation individually
for sec = 1:length(startInds)
    %3a) measure the directionality of the foveations
    
    
    
    secLength = endInds(sec)-startInds(sec);% one less than the total length,
    % but accounted for later, and important for directions between points
    positions = SacStructInput(startInds(sec):endInds(sec),1:2);% the data for this specific foveation
    DX = diff(SacStructInput(startInds(sec):endInds(sec),1)); % displacement at each point for angle calculation
    DY = diff(SacStructInput(startInds(sec):endInds(sec),2));% displacement at each point for angle calculation
        
    %initialize variable for the classification of this section
    secClass = zeros(secLength,1);
    % for each sample, calculate the angle of the movement, and how many
    % windows it is part of
    radians = atan2(DY,DX);
    timesSampled = zeros(secLength,1);
    Rs = zeros(secLength,1);
    %for each section, break it up into 22ms sections with an overlap of
    %6ms, so there is a start at 1, 17, 33 increaseing by 16ms, this
    %translate to sample indices at 1, 9 and 17
    winStarts = [1:8:secLength];
    winEnds = [11:8:secLength, secLength];
    
    for win = 1:length(winEnds)
        timesSampled(winStarts(win):winEnds(win)) =...
            timesSampled(winStarts(win):winEnds(win))+1;
        %calculate the R values
        [Rs(winStarts(win):winEnds(win))] = [Rs(winStarts(win):winEnds(win))+...
            circ_rtest(radians(winStarts(win):winEnds(win)))];
    end
    % to calculate the p value for the rayleigh test, divide the p values 
    % by the number of times the value was calculated (control for overlap) 
    correctedRs = Rs./timesSampled;
    % Now we  separate the sections based on whether the samples are
    % significantly pointed in a direction or not.
    DirSecs = correctedRs<.01; % generates an array with ones when direciton is significant
    
    %check each threshold change, and find sections that are at least 40ms long(20 samples)
    changes = abs(diff(DirSecs));
    sections = [[1; find(changes)-1] [find(changes);length(DirSecs)]];
    %check the length, remove any sections shorter than 20 samples
    goodPortionsInd= sections(diff(sections,1,2)>19,:);
    % get which section number is actually good, so as to acurately assign
    % the threshold crossings
    goodPortionsSec = find(diff(sections,1,2)>19);
    % if there are no sections long enough, but the period is long enough,
    % just calculate for the whole thing
    if isempty(goodPortionsInd)&& length(DirSecs)>20;
        sections = [1 length(DirSecs)];
        goodPortionsInd = [1 length(DirSecs)];
        goodPortionsSec = 1;
    end
    
    %% 4) Calculate Dispersion, Consistency, Path ratio and Spatial Range.
    % if it isn't empty, for each section, compute 4 measures: Dispersion,
    % Consistency, Path ratio and Spatial Range, and test it against the
    % thresholds.
    designation = zeros(size(sections,1),1);
    
    if ~isempty(goodPortionsInd)
        Threshold = zeros(size(sections,1),4);
        for portions = 1:size(goodPortionsInd,1)
            testedPortion = positions(goodPortionsInd(portions,1):goodPortionsInd(portions,2),:);
            % Dispersion 
            
            % First, look at the PCA
            [~, scores] = pca(testedPortion);
            maxPC = max(max(scores(:,1))-min(scores(:,1)),max(scores(:,2))-min(scores(:,2)));
            minPC = min(max(scores(:,1))-min(scores(:,1)),max(scores(:,2))-min(scores(:,2)));
            Dispersion = minPC/maxPC;
            
            %compare the Dispersion to the threshold
            Threshold(goodPortionsSec(portions),1) = Dispersion<DispersionThreshold;
            
            % Consistency
            % calculate the direction consistency, dividing the
            % euclidean distance between the start and end points by
            % the max principal component
            EuclideanDistance = sqrt((testedPortion(end,1)-testedPortion(1,1)).^2+...
                (testedPortion(end,2)-testedPortion(1,2)).^2);
            Consistency = EuclideanDistance/maxPC;
            %compare the consistency to the threshold
            Threshold(goodPortionsSec(portions),2) = Consistency>ConsistDirThreshold;
            
            % Path Ratio
            
            % calculate the ratio between the trajectory length and the
            % Eucliudean Distance.
            Trajectory = sum(sqrt(diff(testedPortion(:,1)).^2+diff(testedPortion(:,2)).^2));
            PathDeviation = EuclideanDistance/Trajectory;
            %compare the PathDeviation to the threshold
            Threshold(goodPortionsSec(portions),3) = PathDeviation>PathRatioThreshold;
            
            %Spatial Range
            
            % calculate the spatial range of the sample
            SpatialRange = sqrt( (max(testedPortion(:,1))-min(testedPortion(:,1))).^2+...
                (max(testedPortion(:,2))-min(testedPortion(:,2) ) ) );
            %compare the spatialRange to the threshold
            Threshold(goodPortionsSec(portions),4) = SpatialRange>SpatialRangeThreshold;
            
            %% 5) Classify the section
            % if all four theresholds are crossed, the section is a Smooth
            % pursuit,
            if sum(Threshold(goodPortionsSec(portions),:))==4
                designation(goodPortionsSec(portions)) = 2;
                % if none of the thresholds are passed, then it is a
                % fixaiton
            elseif sum(Threshold(goodPortionsSec(portions),:))==0
                designation(goodPortionsSec(portions)) = 1;
              
                %if it is undecided, then we have to look at other parts of
                %the section
            end
            
        end
        %% 6) Combine adjacent undetermined sections 
        %check if there were any non-assigned sections
        if any(designation==0)
            % go through and combine any sections that are adjacent and
            % undetermined
            % e.g. if designation = [1 0 0 2 0 1 0 0 0 ], we want to combine 
            % sections corresponding to indices 2 and 3, and 7 8 and 9,
            % with indices 2, 7 and 8 being 'starting' adjacent sections,
            % and 3, 8 and 9 being following sections
            
            % first, find the undesignated sections
            badPors = find(designation==0);
            %get the second part of the adjacent pairs(or more)
            adjacentPors = find(diff(badPors)==1)+1;
            % set the starts of the following adjacent portions to NaN for removal
            sections(badPors(adjacentPors),1) = NaN;
            % set the ends of the starts of the Adjacent portions to NaN
            % for removal
            sections(badPors(adjacentPors)-1,2) = NaN;
            % now, to remove all of the NaNs from Sections, split them up
            % in to the start and end columns, remove the nans, and then
            % recombine them
            secStarts = sections(:,1);
            secStarts(isnan(secStarts)) = [];
            secEnds = sections(:,2);
            secEnds(isnan(secEnds)) = [];
            sections= [secStarts secEnds];
            
            % now, the combined sections are well defined, and all we have
            % to remove is the same indices from the 'designation' array
            designation(badPors(adjacentPors)) = [];
            
            %% 7) Calculate mean direction of all sections, and group with similar sections
            MeanDir = zeros(size(sections,1),1);

            for portions = 1:size(sections,1)
                MeanDir(portions) = circ_mean(radians(sections(portions,1):sections(portions,2)));
            end
            undetermined = find(designation==0);
            %for each undetermined section, find the other sections that
            %have a similar angle (pi/4), and then process them together to
            %calculate a new spatial range
            
            %
            for undetSection = 1:length(undetermined)
                %check the surrounding sections that have similar mean directions
                simDir = find(MeanDir<MeanDir(undetermined(undetSection))+pi/4 & ...
                    MeanDir>MeanDir(undetermined(undetSection))-pi/4 | ...
                    (MeanDir-(2*pi))<MeanDir(undetermined(undetSection))+pi/4 & ...
                    (MeanDir-(2*pi))>MeanDir(undetermined(undetSection))-pi/4 | ...
                    (MeanDir+(2*pi))<MeanDir(undetermined(undetSection))+pi/4 & ...
                    (MeanDir+(2*pi))>MeanDir(undetermined(undetSection))-pi/4);
                testing = [];
                %collect all of the sections that have the same mean dir
                for goodSecs = 1:length(simDir)
                    
                    testing = [testing;positions(sections(simDir(goodSecs),1):sections(...
                        simDir(goodSecs),2),1:2)];
                end
                
                %% 8) Classify based on Path Ratio and Spatial Range
                % now recalculate the Path Ratio of this new section 
                % (the euclidean distance/trajectory)
                positionalDisp = sqrt((testing(end,1)-testing(1,1)).^2+...
                    (testing(end,2)-testing(1,2)).^2)/sum(sqrt(diff(testing(:,1)).^2+diff(testing(:,2)).^2));
                % if the new Path ratio surpasses the threshold, then it is
                % likely a smooth pursuit, and so is compared to the
                % minimum SP spatial range. Otherwise it is more like a
                % fixation, and compared the the max spatial range of a
                % fixation
                if positionalDisp>PathRatioThreshold
                      
                    % if it looks like an SP, then check the spatial range
                    % of this section with the whole testing bit
                    NewSpatialRange = sqrt( (max(testing(:,1))-min(testing(:,1))).^2+...
                        (max(testing(:,2))-min(testing(:,2) ) ) );
                    % if the spatial range is above the minimum range for a smooth
                    % pursuit, then it is part of a larger smooth pursuit,
                    % otherwise, it is a fixation.
                    if NewSpatialRange>MinSpatialRange
                        designation(undetermined(undetSection)) = 2;
                    else
                        designation(undetermined(undetSection)) = 1;
                    end
                else
                    %if it isn't like a smooth pursuit in terms of the
                    %pathRatioThreshold, then test it with the spatial
                    %range, see if it is too big to be a fixation,
                    %otherwise it is a fixation
                    SpatialRange = sqrt( (max(testing(:,1))-min(testing(:,1))).^2+...
                        (max(testing(:,2))-min(testing(:,2) ) ) );
                    %test if the new spatial range is too great to be a
                    %fixation
                    if SpatialRange > SpatialRangeThreshold
                        designation(undetermined(undetSection)) = 2;
                    else % otherwise it is a fixation
                        designation(undetermined(undetSection)) = 1;
                    end
                end
            end
            
        end
        for finalSect = 1:size(sections,1)
            secClass(sections(finalSect,1):sections(finalSect,2))=...
                designation(finalSect);
        end

    end
    secClass = [secClass;secClass(end)];
    classification(startInds(sec):endInds(sec)) = secClass;
end
% figure
% 
% scatter(SacStructInput(classification==2,1),SacStructInput(classification==2,2),'b')
% hold on;scatter(SacStructInput(classification==1,1),SacStructInput(classification==1,2),'r')
% scatter(SacStructInput(classification==0,1),SacStructInput(classification==0,2),'g')
% axis([-20 20 -15 15]);
% 
% figure
% hold on
% scatter(SacStructInput(classification==2,3),SacStructInput(classification==2,1),'b')
% scatter(SacStructInput(classification==1,3),SacStructInput(classification==1,1),'r')
% scatter(SacStructInput(classification==0,3),SacStructInput(classification==0,1),'g')
% plot(SacStructInput(:,3),SacStructInput(:,1))
% scatter(SacStructInput(classification==2,3),SacStructInput(classification==2,2),'b')
% scatter(SacStructInput(classification==1,3),SacStructInput(classification==1,2),'r')
% scatter(SacStructInput(classification==0,3),SacStructInput(classification==0,2),'g')
% plot(SacStructInput(:,3),SacStructInput(:,2))
% 
% boop=1;
% boop=1;
end


