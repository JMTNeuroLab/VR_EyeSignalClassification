
function [startPoint, endPoint,peakVel,peakPoint] = calcEndPoints(pos)
%% 
% this script takes the potential saccadic periods found by
% SaccadeProcessing and finds the start point, end point, peak velocity and
% the time of peak velocity (peakPoint).
% This is modified for primates from Larsson et al. (2013)
% 
% this processes each saccade independently
% 1) Calculate the Main direction and the changes in direction
% 2) 

%% initialize variables
p=0; % change to 1 if figures are desired
MainDirThresh = 20;
SampleChangeThresh = 20;
startPoint = NaN;
endPoint = NaN;
peakVel = NaN;
peakPoint = NaN;

saccade = 1;
while saccade == 1
    %% 1) find the main direction
    % calculate the change in position at each point
    DX = diff(pos(:,1));
    
    DY = diff(pos(:,2));
    % use the change in position to generate an angle at each point
    Angle = atan2d(DY,DX)+180;
    % calculat ethe velocity at each point
    vel = sqrt(diff(pos(:,1)).^2+diff(pos(:,2)).^2)/.002;
    
    [peakVel,peak] = max(vel);
    if peakVel<50
        saccade = 0;
        continue
    end
    % the velThresh is the greater of 30deg/s or 1/5 of the peakVel
    velThresh = max(peakVel*.2,30);
    
    if peak<8||peak>length(pos)-20 % check that there is enough space to 
        % calculate the endpoints, if not, skip this saccade
        saccade = 2;
        continue
    end
    % length(AngleStart)+length(AngleEnd)
    % length(pos)
    % the main direction, which is the mean over three points at the peak
    MainDir = meanangle(Angle(peak-1:peak+1));   
    
    % calculate the difference between the direction at a point and the
    % MainDir for the angles before the peak and after the peak
    AngleStart = min(mod(Angle(1:peak)-MainDir, 360),mod(MainDir - Angle(1:peak), 360));
    AngleEnd = min(mod(Angle(peak:end)-MainDir,360),mod(MainDir - Angle(peak:end), 360));
    %% plots
    if p==1
        figure('position',[200 200 1500 800])
        subplot(2,3,1)
        plot(pos(2:end,1),pos(2:end,2),'b-',pos(2:end,1),pos(2:end,2),'bx','LineWidth',3)
        xlabel('Horizontal Eye Position X (deg)')
        ylabel('Vertical Eye Position Y (deg)')
        title('Eye Position on Screen')
        ax = gca;
            ax.FontSize = 16;

        subplot(2,3,2)
        
        
        plot(pos(2:end,1),'LineWidth',3)
        xlabel('Samples (2ms)')
        ylabel('X Position (deg)')
        title('Eye Position in X Axis Over Time')
        ax = gca;
            ax.FontSize = 16;
    end
    
    %% 2) Find the MainDir threshold crossings for onset and offset
    %see if there are more than three samples (6ms) where the direction is more than
    %20 deg different than the main direction, or if it is 60 deg different
    %for one samples. Find the last time that this happens (closest to the
    %peak)
    
    %initialize default values
    MainDirStart = NaN;
    MainDirEnd = NaN;
    if any(AngleStart>MainDirThresh)
        %calculate startpoint, either above the threshold for 6 ms, or
        %above 3 times the threshold once
        singleMainDirStart = max([find(AngleStart>(3*MainDirThresh),1,'last'),NaN]);
        tripleMainDirStart = max([peak-strfind([flipud(AngleStart>MainDirThresh)]',[1,1,1]),NaN]);
        MainDirStart = max(tripleMainDirStart,singleMainDirStart);
        % check that the velocity is lower than thevelocity threshold,
        % otherwise, set it at the first point below the velocity threshold
         if ~isnan(MainDirStart) && vel(MainDirStart)>velThresh
            MainDirStart = find(vel(1:peak)<velThresh,1,'last');
            if isempty(MainDirStart)
                saccade = 0;
                continue
            end
        end
        if p==1 && ~isnan(MainDirStart) %if plotting, plot this line
            line([MainDirStart MainDirStart],[min(pos(:,1)-2),max(pos(:,1)+2)],'color','r','LineWidth',3)
        end
    end
    %repeat for the end point, however, use the first time this happens
    if any(AngleEnd>MainDirThresh)
        %calculate endpoint
        singleMainDirEnd = min([find(AngleEnd>3*MainDirThresh,1),NaN]);
        tripleMainDirEnd = min([strfind([AngleEnd>MainDirThresh]',[1 1 1]),NaN]);
        MainDirEnd = min(singleMainDirEnd,tripleMainDirEnd);
         % if the end point is the last point, then we add the
         % final velocity to the end for the velocity check
        if ~isnan(MainDirEnd) && MainDirEnd+peak>length(vel)
            vel = [vel;vel(end)];
        end
        % check that the velocity is lower than thevelocity threshold,
        % otherwise, set it at the first point below the velocity threshold
        if ~isnan(MainDirEnd)&& vel(MainDirEnd+peak)>velThresh 
            MainDirEnd = find(vel(peak:end)<velThresh,1);
            if isempty(MainDirEnd)
                saccade = 0;
                continue
            end
        end
        if p==1 && ~isnan(MainDirEnd)%if plotting, plot this line
            line([MainDirEnd+peak MainDirEnd+peak],[min(pos(:,1)-2),max(pos(:,1)+2)],'color','r','LineWidth',3)
        end
    end
    %% 3) Find the sample-to-sample direction threshold crossings
    % get the smallest angle between each point
    SampleDirectionStart = min(abs(diff(Angle(1:peak))),...
        min(abs(diff(Angle(1:peak))-360),abs(diff(Angle(1:peak))+360)));
    SampleDirectionEnd = min(abs(diff(Angle(peak:end))),...
        min(abs(diff(Angle(peak:end))-360),abs(diff(Angle(peak:end))+360)));
    
    %initialize default values
    IncDirStart = NaN;
    IncDirEnd = NaN;
    
     % see if there are more than three samples (6ms) where the direction is
     % always changing by more than 20 deg, or if there is a single change of
     % more than 60 deg. 
     % Find the first time that this happens (closest to the peak) for the
     % start
    if any(SampleDirectionStart>SampleChangeThresh)
        %
        singleIncDirStart = max([find(SampleDirectionStart>3*SampleChangeThresh,1,'last'),NaN]);
        tripleIncDirStart = max([peak-strfind([flipud(SampleDirectionStart>...
            SampleChangeThresh)]',[1 1 1]),NaN]);
        IncDirStart = max(singleIncDirStart,tripleIncDirStart);
        % check that the velocity is lower than thevelocity threshold,
        % otherwise, set it at the first point below the velocity threshold
        if ~isnan(IncDirStart)&& vel(IncDirStart)>velThresh 
            IncDirStart = find(vel(1:peak)<velThresh,1,'last');
            if isempty(IncDirStart)
                saccade = 0;
                continue
            end
        end
        if p==1 &&~isnan(IncDirStart)%if plotting, plot this line
            line([IncDirStart IncDirStart],[min(pos(:,1)-2),max(pos(:,1)+2)],'color','g','LineWidth',3)
        end
    end
    % do the same to find the endpoint
    if any(SampleDirectionEnd>SampleChangeThresh)
        singleIncDirEnd = max([find(SampleDirectionEnd>SampleChangeThresh*3,1),NaN]);
        tripleIncDirEnd = min([strfind([SampleDirectionEnd>SampleChangeThresh]',[1 1 1]),NaN]);
        IncDirEnd = min(singleIncDirEnd,tripleIncDirEnd);
        % if the end point is the last point, then we add the
         % final velocity to the end for the velocity check
        if ~isnan(IncDirEnd) && IncDirEnd+peak>length(vel)
            vel = [vel;vel(end)];
        end
        % check that the velocity is lower than thevelocity threshold,
        % otherwise, set it at the first point below the velocity threshold
        if ~isnan(IncDirEnd) && vel(IncDirEnd+peak)>velThresh
            IncDirEnd = find(vel(peak:end)<velThresh,1);
            if isempty(IncDirEnd)
                saccade = 0;
                continue
            end
        end
        if p==1 && ~isnan(IncDirEnd)%if plotting, plot this line
            line([IncDirEnd+peak IncDirEnd+peak],[min(pos(:,1)-2),max(pos(:,1)+2)],'color','g','LineWidth',3)
        end
        
    end
    %% plotting
    if p==1
        %% plot Y
        subplot(2,3,3)
        plot(pos(2:end,2),'LineWidth',3)
        xlabel('Samples (2ms)')
        ylabel('Position in Y Axis (deg)')
        title('Eye Position in Y Axis Over Time')
        if any(AngleStart>MainDirThresh)
            %calculate startpoint
            line([MainDirStart MainDirStart],[min(pos(:,2)-2),max(pos(:,2)+2)],'color','r','LineWidth',3)
        end
        if any(AngleEnd>MainDirThresh)
            %calculate endpoint
            line([MainDirEnd+peak MainDirEnd+peak],[min(pos(:,2)-2),max(pos(:,2)+2)],'color','r','LineWidth',3)
        end
        if any(SampleDirectionStart>SampleChangeThresh);
            line([IncDirStart IncDirStart],[min(pos(:,2)-2),max(pos(:,2)+2)],'color','g','LineWidth',3)
        end
        if any(SampleDirectionEnd>SampleChangeThresh)
            line([IncDirEnd+peak IncDirEnd+peak],[min(pos(:,2)-2),max(pos(:,2)+2)],'color','g','LineWidth',3)
        end
        ax = gca;
            ax.FontSize = 16;
        %% plot Velocity
        subplot(2,3,4)
        plot(vel,'LineWidth',3)
        xlabel('Samples (2ms)')
        ylabel('Angular Velocity of Eye Position (deg/s)')
        title('Angular Velocity of Eye Position Over Time')
        if any(AngleStart>MainDirThresh)
            line([MainDirStart MainDirStart],[min(sqrt(diff(pos(:,1)).^2+diff(pos(:,2)).^2)/.002-2),...
                max(sqrt(diff(pos(:,1)).^2+diff(pos(:,2)).^2)/.002+2)],'color','r','LineWidth',3)
        end
        if any(AngleEnd>MainDirThresh)
            line([MainDirEnd+peak MainDirEnd+peak],[min(sqrt(diff(pos(:,1)).^2+diff(pos(:,2)).^2)/.002-2),...
                max(sqrt(diff(pos(:,1)).^2+diff(pos(:,2)).^2)/.002+2)],'color','r','LineWidth',3)
        end
        if any(SampleDirectionStart>SampleChangeThresh);
            line([IncDirStart IncDirStart],[min(sqrt(diff(pos(:,1)).^2+diff(pos(:,2)).^2)/.002-2),...
                max(sqrt(diff(pos(:,1)).^2+diff(pos(:,2)).^2)/.002+2)],'color','g','LineWidth',3)
        end
        if any(SampleDirectionEnd>SampleChangeThresh)
            line([IncDirEnd+peak IncDirEnd+peak],[min(sqrt(diff(pos(:,1)).^2+diff(pos(:,2)).^2)/.002-2),...
                max(sqrt(diff(pos(:,1)).^2+diff(pos(:,2)).^2)/.002+2)],'color','g','LineWidth',3)
        end
        ax = gca;
             ax.FontSize = 16;
         ax.YLim(1) = 0;
        subplot(2,3,5)
        plot([AngleStart;AngleEnd],'r','LineWidth',3)
        xlabel('Samples (2ms)')
        ylabel({'Absolute Difference from','the Main Direction (deg)'})
        title({'Absolute Difference from the Main direction','of the Saccade Over Time'})
        ax = gca;
            ax.FontSize = 16;
            
        subplot(2,3,6)
        plot([SampleDirectionStart;SampleDirectionEnd],'g','LineWidth',3)
        xlabel('Samples (2ms)')
        ylabel({'Absolute Sample to Sample Change',' in Direction (deg)'})
        title({'Absolute Sample to Sample Change',' in Direction Over Time'})
        ax = gca;
            ax.FontSize = 16;
        subplot(2,3,1)
        
        hold on
        %draw a line for scale that is 2 degrees to the left of the saccade trace,
        %and 6 degrees along the y axis
        line([min(pos(:,1))-2 min(pos(:,1))-2],...
            [mean([max(pos(:,2)),min(pos(:,2))])-3 mean([max(pos(:,2)),min(pos(:,2))])+3],'LineWidth',3);
        colors = [0 1 0;1 0 0];
        scatter([pos(1,1),pos(end,1)],...
            [pos(1,2),pos(end,2)],100,colors,'d','LineWidth',3)
        if ~isnan(IncDirStart)&& ~isnan(IncDirEnd)
            scatter([pos(IncDirStart,1),pos(peak+IncDirEnd,1)],...
                [pos(IncDirStart,2),pos(peak+IncDirEnd,2)],10,'g','LineWidth',3)
        end
        if ~isnan(MainDirStart) && ~isnan(MainDirEnd)
            scatter([pos(MainDirStart,1),pos(peak+MainDirEnd,1)],...
                [pos(MainDirStart,2),pos(peak+MainDirEnd,2)],10,'r','LineWidth',3)
        end
         pause
        close
    end%other Plots, cancelled if p==0
    %% 5) Define the saccade onset and offset 
    % Check which of the direction thresholds is closest to the peak for 
    % both the start and the end points. 
    if isnan(MainDirStart) && isnan(IncDirStart)
        bad =1;
    else
        % add one at the start because the threshold crossing identifies
        % the direction change, which looses a point during direction
        % calculation with diff
        MainDirStart = MainDirStart+1;
        IncDirStart = IncDirStart+1;
        startPoint = pos(max(MainDirStart,IncDirStart),3);
    end
    peakPoint = pos(peak,3);
    
    if isnan(MainDirEnd) && isnan(IncDirEnd)
        bad = 2;
    else
        % add the peak +1 to the threshold because the search starts at the
        % peak, and the direction search looses a point during direction
        % calculation with diff
        MainDirEnd = min(MainDirEnd+1+peak,length(pos));
        IncDirEnd = min(IncDirEnd+1+peak,length(pos));
        endPoint = pos(min(MainDirEnd,IncDirEnd),3);
    end
    % check that saccade was at least 10ms
    if endPoint-startPoint<.01
        bad = 3;
    end
    saccade = 2;
end
end

%%



