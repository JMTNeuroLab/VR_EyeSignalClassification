function [SacStructOut] = HcTask_SaccadeProcessing(SacStructInput,PreSmoothed,toPlot)
% This script is used to classify eye movement data collected at 500hz in a
% macaque viewing a dynamic display.
% A matrix with the eye-position data in degrees visual angle
% is input, with columns (x, y, time(s)) and the output is a 
% structure with eye movement classification as well as saccade kinematics.
% The 'PreSmoothed' variable allows you to process data that has already
% been smoothed, or unsmoothed data. toPlot allows for the trial data to be
% plotted with the overlayed classification.

% These algorithms are based on Larsson et al. (2013 and 2015), and are
% used in Corrigan et al (2017) and Doucet et al, (2017).
% Also used are functions from Philip Berens' circular statistics toolbox
% and Jeff Dunne (meanaAngle).

% The basic stucture of this script is that it 
% 1) finds an acceleration threshold for the x and y dimensions
% 2) find periods that are potential saccades
% 3) Process each potential saccadic period (calcEndPoints.m)
% 4) Gather Saccade kinematics info
% 5) Test for Post-Saccadic Oscillations (PSO's)
% 6) generate an index of saccade and PSO classification
% 7) classify the foveations (FixVsSPAnalysis.m)
% 8) plot the figures if desired

%to plot the example, go to calcEndPoints.m, and change p to 1
%   SacStructInput =BigStruct.(Types{sessionType}).Trials.(trialIDs{trl}).PositionMatrix_degOnScreen(:,4:6)
%   SacStructInput = XMazeStruct.Trials.
%   SacStructInput=BigStruct.(Types{sessionType}).Trials.(trialIDs{trl}).EyeDegrees
%   SacStructInput=BigStruct.(Types{sessionType}).Trials.(trialIDs{trl}).NonSmoothedEyes

%% initialize variables
%number of samples to exceed to count as a saccade 
durationThresh = 4;% 8ms, so saccades must be 10ms

excessEstimate = 1000;
%the cells to hold the processed data before putting it into the array
TotalStartTime = NaN(excessEstimate,1);
TotalEndTime = NaN(excessEstimate,1);
TotalPeakTime = NaN(excessEstimate,1);
TotalPeakVelocity = NaN(excessEstimate,1);
StartPointX  = NaN(excessEstimate,1);
StartPointY = NaN(excessEstimate,1);
EndPointX  = NaN(excessEstimate,1);
EndPointY  = NaN(excessEstimate,1);
Amplitude  = NaN(excessEstimate,1);
Duration  = NaN(excessEstimate,1);
Direction  = NaN(excessEstimate,1);
PSOEndTotal = NaN(excessEstimate,1);
%the window for smoothing
wind = 20*[ -1 -1 -1 -1 0 1 1 1 1]; % window is dependent on sampling rate,
% and adapted from Larsson et al. (2013) to calculate velocity in
% degrees/sec
% Tsac keeps track of the total saccades
Tsac = 0;
saccadesExist = 1;
red = false(length(SacStructInput),1); %fix inds logical
blue = false(length(SacStructInput),1); %smooth pursuit logical
green =  false(length(SacStructInput),1); %saccades logical
magenta = false(length(SacStructInput),1);
EyeDataCopy = SacStructInput;
%%
while saccadesExist ==1
    %% 1) Calculate the velocity Threshold
    %the position has already been smoothed, so just calculate the velocity
    %by the difference of the two vectors
    
    %Calculate the Velocity and the Acceleration for the trial
    if PreSmoothed ==1
        Vel = abs(diff(SacStructInput(:,1:2))./...
            [diff(SacStructInput (:,3)),diff(SacStructInput(:,3))]);
    else
        Vel = abs([conv(SacStructInput(:,1),wind,'same'),conv(SacStructInput(:,2),wind,'same')]);
    end
    Acc = abs([conv(Vel(:,1),wind,'same'),conv(Vel(:,2),wind,'same')]);
    
    %Calculate the threshold iteratively, taken from Larsson, L., Nystrom,
    %     M., & Stridh, M. (2013). Detection of saccades and postsaccadic
    %     oscillations in the presence of smooth pursuit. IEEE Transactions on
    %     Biomedical Engineering, 60(9), 2484?2493.
    %     http://doi.org/10.1109/TBME.2013.2258918
    
    %   goes through and gets the mean of all the points below the threshold,
    %   calculates a new threshold based on the new mean + 6* the standard
    %   deviation, and then iterates again. stops after it stops decreasing
    %   by more than 3 deg/s^2
    
    %initialize the values for the differences between iterations of the
    %threshold calculator, and an initial Peak Threshold value
    Threshdiff = [3,3];
    PeakThresh = [10000,10000];
    for var = 1:2
        
        while Threshdiff(var)>1
            MeanVel(var) = mean(Acc(Acc(:,var)<PeakThresh(var),var));
            STD(var) = std(Acc(Acc(:,var)<PeakThresh(var),var));
            NewThresh = MeanVel(var)+(6*STD(var));
            Threshdiff(var) = PeakThresh(var)-NewThresh;
            PeakThresh(var) = NewThresh;
            
        end
        
    end
    %% 2) find periods that are potential saccades
    %define the portions that are above the threshold for each component
    isSac = Acc(:,1)>PeakThresh(1) | Acc(:,2)>PeakThresh(2);
    if isempty(isSac) %if there are none, exit, as there is nothing to analyze
        saccadesExist = 0;
        continue
    end
    if diff(isSac(1:2))==-1 || isSac(1)==1
        %NonSacStart is all the points where isSac changes from 1 to 0, but
        %we have to check that it doesn't do this at the first index
        NonSacStart = find(diff(isSac)==-1);
    else     NonSacStart = [1;find(diff(isSac)==-1)];
    end
    %NonSacEnd is at the points where a sacade starts, so where isSac
    %starts to go to 1
    NonSacEnd = [find(diff(isSac)==1);length(isSac)];
    %make a matrix and find the length of the subthreshold periods
    if isempty(NonSacStart)||isempty(NonSacEnd)
        saccadesExist = 0;
        continue
    end
    startsandEndsTemp = [NonSacStart NonSacEnd(1:length(NonSacStart))];
    startsandEndsTemp(:,3) = diff(startsandEndsTemp,1,2);
    %if the span of sub threshold is too short, remove the start and end points
    tooShort = find(startsandEndsTemp(:,3)<20);
    NonSacStart(tooShort) = [];
    NonSacEnd(tooShort) = [];
    %now to hiccups, where there isn't actually a saccade, because
    %it is too short
    if isempty(NonSacStart(2:end))||isempty(NonSacEnd(1:length(NonSacStart)-1))
        saccadesExist = 0;
        continue
    end
    startsandEnds = [NonSacEnd(1:length(NonSacStart)-1) NonSacStart(2:end)];
    startsandEnds(:,3) = diff(startsandEnds,1,2);
    %get the times that are longer than 10ms
    SaccadeTimes = find(startsandEnds(:,3)>durationThresh);
    if isempty(SaccadeTimes)
        saccadesExist = 0;
        continue
    end
    %use these as the start and end of the bad data
    %check that there is some room for the padding, remove saccades that
    %start too early, they canont be properly analyzed
    if startsandEnds(SaccadeTimes(1),1)<5
        %check that the saccade doesn't end late and have to be removed as
        %well
        if startsandEnds(SaccadeTimes(end),2)+20>length(SacStructInput)
            SacEnd = [startsandEnds(SaccadeTimes(2:end-1),2)]+20;
            SacStart = [startsandEnds(SaccadeTimes(2:end-1),1)]-4;
        else SacEnd = startsandEnds(SaccadeTimes(2:end),2)+20;
            SacStart = [startsandEnds(SaccadeTimes(2:end),1)]-4;
        end
        %if the start doesn't need to be modified, still check the end
    else
        if startsandEnds(SaccadeTimes(end),2)+20>length(SacStructInput)
            SacEnd = [startsandEnds(SaccadeTimes(1:end-1),2)]+20;
            SacStart = startsandEnds(SaccadeTimes(1:end-1),1)-4;
            
        else SacEnd = startsandEnds(SaccadeTimes,2)+20;
            SacStart = startsandEnds(SaccadeTimes,1)-4;
        end
    end
    %% 3) Process each potential saccadic period
    %initialize variables
    StartPoint = NaN(length(SacStart),1);
    EndPoint = NaN(length(SacStart),1);
    Peak = NaN(length(SacStart),1);
    PeakTime =  NaN(length(SacStart),1);
    %actual saccade counter
    saccadeCount = 0;
    %go through all of the  possible saccade periods
    for stamp = 1:length(SacStart)
        
        %make sure there aren't any NaNs in the sample, and if not,
        %increase the saccade counter by one, and calculate the pertinent
        %Points
        if ~any(isnan(SacStructInput(SacStart(stamp):SacEnd(stamp))))
            saccadeCount = saccadeCount+1;
            [StartPoint(saccadeCount),EndPoint(saccadeCount),Peak(saccadeCount),...
                PeakTime(saccadeCount)] = calcEndPoints(SacStructInput...
                (SacStart(stamp):SacEnd(stamp),:));
        end
    end
    %remove any NaN from bad eye samples
    toremove = (isnan(StartPoint)|isnan(EndPoint)|isnan(Peak)); % find any saccades that have NaN in them
%     RemovedStart = StartPoint(toremove); % these can be used to analyse
%     how many and when saccades are removed
%     RemovedEndPoint = EndPoint(toremove);
    StartPoint(toremove) = []; % remove this info from the indicies
    EndPoint(toremove) = [];
    Peak(toremove) = [];
    PeakTime(toremove) = [];
    PSOEndTotal = NaN(length(StartPoint),1);
    %% 4) Gather Saccade kinematics info
    % Cycle through Saccades
    for sac = 1:length(StartPoint)
        Tsac = Tsac+1; % increase the total saccade count
        % the following adds the data from the current saccade to the
        % arrays for different data types
        TotalStartTime(Tsac) = StartPoint(sac); % ad
        TotalEndTime(Tsac) = EndPoint(sac);
        TotalPeakTime(Tsac) = PeakTime(sac);
        TotalPeakVelocity(Tsac) = Peak(sac);
        % use inddexing to get the position of the start/end of the saccade
        StartPointX(Tsac) = SacStructInput(SacStructInput(:,3)==StartPoint(sac),1); % xposition of the start of the saccade
        StartPointY(Tsac) = SacStructInput(SacStructInput(:,3)==StartPoint(sac),2);
        EndPointX(Tsac) = SacStructInput(SacStructInput(:,3)==EndPoint(sac),1);
        EndPointY(Tsac) = SacStructInput(SacStructInput(:,3)==EndPoint(sac),2);
        Amplitude(Tsac) = sqrt((StartPointX(Tsac)-EndPointX(Tsac)).^2 + ...
            (StartPointY(Tsac)-EndPointY(Tsac)).^2);
        Duration(Tsac) = EndPoint(sac)-StartPoint(sac);
        Direction(Tsac) = atan2d((EndPointY(Tsac)-StartPointY(Tsac)),...
            (EndPointX(Tsac)-StartPointX(Tsac)));
        %% 5) Test for Post-Saccadic Oscillations (PSO's)
        
        % if the saccade is the last saccade, there might not be enough
        % time left to use the initial 80 ms used, but this isn't
        % necessary. If there are less then 38 samples (78ms) then use the
        % end of the dataset. Just check that there are more than 12
        % samples before proceeding
        %CheckX
        if sac+1>length(StartPoint)
            calcLength = min(find(SacStructInput(:,3)==EndPoint(sac))+38,length(SacStructInput));
        else
            calcLength = min([find(SacStructInput(:,3)==EndPoint(sac))+38,...
                find(SacStructInput(:,3)==StartPoint(sac+1))-1]);
        end
        if calcLength >10 && [calcLength - find(SacStructInput(:,3)==EndPoint(sac))]>10
            PPSOX = SacStructInput(find(SacStructInput(:,3)==EndPoint(sac)):...
                calcLength,1:2:3);
            
            [PSOEndX,Rend] = PSOCalc(PPSOX);
            %CheckY
            PPSOY = SacStructInput(find(SacStructInput(:,3)==EndPoint(sac)):...
                calcLength,2:3);
            
            [PSOEndY, REnd] = PSOCalc(PPSOY);
            % use the longer of the two PSO durations
            PSOEndTotal(sac) = max(PSOEndX,PSOEndY);
            
        end
        %% 6) generate an index of saccade and PSO classification
        %define points where the signal is a saccade
        green(SacStructInput(:,3)>=StartPoint(sac)&SacStructInput(:,3)<=EndPoint(sac)) = 1;
        %define points where the signal is a PSO. Add on sample (2ms) to the start to
        %ensure that a point that is considered a saccade is not included
        %in the PSO
        magenta(SacStructInput(:,3)>=(EndPoint(sac)+.002)&SacStructInput(:,3)<=PSOEndTotal(sac)) = 1;
        % nan all the defined points so that only the foveations are left
        % for further processing
        SacStructInput(SacStructInput(:,3)>=StartPoint(sac)&...
            SacStructInput(:,3)<=PSOEndTotal(sac),1:2) = NaN;
        %     figure
    end
    
    saccadesExist = 2;
end
% remove any NaN data from the data arrays
TotalStartTime(isnan(TotalStartTime)) = [];
TotalEndTime(isnan(TotalEndTime)) = [];
TotalPeakTime(isnan(TotalPeakTime)) = [];
TotalPeakVelocity(isnan(TotalPeakVelocity)) = [];
StartPointX(isnan(StartPointX)) = [];
StartPointY(isnan(StartPointY)) = [];
EndPointX(isnan(EndPointX)) = [];
EndPointY(isnan(EndPointY)) = [];
Amplitude(isnan(Amplitude)) = [];
Duration(isnan(Duration)) = [];
Direction(isnan(Direction)) = [];
PSOEndTotal(isnan(PSOEndTotal)) = [];

%% 7) classify the foveations
classifications = FixVsSPAnalysis(SacStructInput);
blue(classifications==2) = true;
red(classifications==1) = true;

magenta = logical(magenta);
green = logical(green);
black = isnan(EyeDataCopy(:,1));

  %% 8) plot figures
  if exist('toPlot','var') && toPlot ==1
      figure
      subplot(2,1,1)
      hold on
      bigdata = [EyeDataCopy(green,3),EyeDataCopy(green,1),ones(length(EyeDataCopy(green)),1);...
          EyeDataCopy(magenta,3),EyeDataCopy(magenta,1),ones(length(EyeDataCopy(magenta)),1).*4;...
          EyeDataCopy(red,3),EyeDataCopy(red,1),ones(length(EyeDataCopy(red,1)),1).*2;...
          EyeDataCopy(blue,3),EyeDataCopy(blue,1),ones(length(EyeDataCopy(blue)),1).*3; ];
      gscatter(bigdata(:,1),bigdata(:,2),bigdata(:,3),'grbm','oooo',[5,5,5,5])
      plot(EyeDataCopy(:,3),EyeDataCopy(:,1),'k','LineWidth',1)
      
      %
      legend({'Saccade','Fixation', 'Smooth Pusuit','Post-saccadic Osscilation'    })
      ylabel('Degrees visual angle X')
      bigfig = gca;
      trialLength = EyeDataCopy(end)-EyeDataCopy(1,3);
      quarters = round(length(EyeDataCopy(:,3))/4);
      bigfig.XTick = [EyeDataCopy(1,3) EyeDataCopy(quarters,3) ...
          EyeDataCopy((2*quarters),3) EyeDataCopy(3*quarters,3) EyeDataCopy(end,3)];
      bigfig.XTickLabel= {'0', num2str(length(EyeDataCopy(:,1))/2) num2str(length(EyeDataCopy(:,1))) ...
          num2str(length(EyeDataCopy(:,1))*1.5) num2str(length(EyeDataCopy(:,1))*2)};
      bigfig.YLim = [-20 20];
      xlims = bigfig.XLim;
      subplot(2,1,2)
      hold on
      scatter(EyeDataCopy(green,3),EyeDataCopy(green,2),20,'go','LineWidth',2)
      scatter(EyeDataCopy(red,3),EyeDataCopy(red,2),20,'ro')
      scatter(EyeDataCopy(blue,3),EyeDataCopy(blue,2),20,'bo')
      scatter(EyeDataCopy(magenta,3),EyeDataCopy(magenta,2),20,'mo')
      plot(EyeDataCopy(:,3),EyeDataCopy(:,2),'k','LineWidth',1)
      % legend({'Fixation', 'Smooth Pusuit','Post-saccadic Osscilation',...
      %             'Saccade'})
      xlabel('Time (ms)')
      ylabel('Degrees visual angle Y')
      bigfig = gca;
      trialLength = EyeDataCopy(end)-EyeDataCopy(1,3);
      quarters = round(length(EyeDataCopy(:,3))/4);
      bigfig.XTick = [EyeDataCopy(1,3) EyeDataCopy(quarters,3) ...
          EyeDataCopy((2*quarters),3) EyeDataCopy(3*quarters,3) EyeDataCopy(end,3)];
      bigfig.XTickLabel= {'0', num2str(length(EyeDataCopy(:,1))/2) num2str(length(EyeDataCopy(:,1))) ...
          num2str(length(EyeDataCopy(:,1))*1.5) num2str(length(EyeDataCopy(:,1))*2)};
      bigfig.YLim = [-15 15];
      bigfig.XLim = xlims;
  end

%% Output
SacStructOut = struct('StartTime',TotalStartTime,'EndTime',TotalEndTime,...
    'PeakTime',TotalPeakTime,'PeakVelocity',TotalPeakVelocity,'Duration',...
    Duration,'Amplitude',Amplitude,'StartPointX',StartPointX,'EndPointX',...
    EndPointX,'StartPointY',StartPointY,'EndPointY',EndPointY,...
    'Direction',Direction,'PostSaccadicOscillationEnd',PSOEndTotal,'Saccade'...
    ,green,'PSO',magenta,'Fixation',red,'SmoothPursuit',blue,'OffScreen',black);
end