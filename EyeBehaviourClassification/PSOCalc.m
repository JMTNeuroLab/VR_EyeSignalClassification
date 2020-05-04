function [PSOendpoint,REnd] = PSOCalc(Input)
%takes a section of the post saccaddic period( X or Y) and finds if there
%is an oscillation in it, and returns the timestamp of the end of the oscillation.
% 
% 1) Find where the linear part of the signal starts
% 2) F


PPSO = Input(:,1);
Time = Input(:,2);

%%
%initialize the variables
isPSO = 1; %changed if the section is not a PSO, or to end the algorithm
ThreshAngle = .08; %the minimum difference between the modeled signal and the decay
ender = min(20,length(PPSO)-7); %intially, set endpoint at 40 ms, or 20 samples, may be extended later
GNended = 0; %the variable for the while loop that checks if an end is found for GN
Vender = ender; %the variable ender that moves backwards until a linear signal is found in GN
PSOendpoint = Time(1);
REnd = 'good';
removed = 0; %logs any changes in the beginning of the signal (such as cutting it to make the model fit)
%check that the signal has evened out - so that the slope before and after
%have the same sign
while isPSO==1
    %% 1) Find where the linear part of the signal starts
    TLrefslope = polyfit([1:5]'*.002,PPSO(ender:ender+4),1);
    if ender-7>0 %check that ender is long enough to have a slope, if it is
        % too short, there must be not enough samples, so don't check the
        % slope, just keep ender there
        TLtestslope = polyfit([1:5]'*.002,PPSO(ender-4:ender),1);
        %check if there is any NaN in the sample
        if any(isnan(PPSO(1:ender)))
            isPSO = 0;
            continue
        end
        
        %check if we have to extend the sample (if there are opposite signs of
        %the slopes of the lines)
        if TLtestslope(1)*TLrefslope(1)<0
            ender = min(ender+10,length(PPSO));
        end
    end
    %check if there is any NaN in the sample
    if any(isnan(PPSO(1:ender)))
        isPSO = 0;
        continue
    end
    %Go through and find out where the linear part starts
    Vender = ender;
    while GNended == 0
        GNrefslope = polyfit([1:ender-Vender+3]'*.002, PPSO(Vender-2:ender),1)';
        GNtestslope = polyfit([1:2]'*.002, PPSO(Vender-3:Vender-2),1)';
%       slopes = [slopes GNrefslope(1)-GNtestslope(1)];
        if abs(GNrefslope(1)-GNtestslope(1))>17
            GNended = 1;
        else
            Vender = Vender-1;
        end
        if Vender<=5
            GNended = 1;
            isPSO = 0;
            continue
        end
        
    end
    %% 2) Transpose the signal so that the start of the linear part is 0, and
    % set the signal to the right of this to 0 as well
    GN = ones(ender,1)*PPSO(Vender-1);
    GN(1:Vender-1) = PPSO(1:Vender-1);
    GN = GN-PPSO(Vender-1);
    % check if there is anything left of the signal
    if isempty(find(GN,1))
        isPSO = 0;
        continue
    end
    %Check if there is enough of an amplitude to warrant a PSO
    if max(abs(GN))<.15
        isPSO = 0;
        continue
    end
    %% 3) Calculate the Impulse models for 1:4 poles
       
    adequateRMSE = 0; % variable to stop while loop
    
    while adequateRMSE ==0;
        GI = zeros(length(GN),4); % initialize the matrix to hold the impulse decay patterns
        RMSE = [1,1,1,1];
        for test = 1:4
            [num{test},denom{test}] = prony(GN,0,test);
            GI(:,test) = impz(num{test},denom{test},length(GN));
            RMSE(test) = (sqrt(sum((GN(:)-GI(:,test)).^2)/numel(GN)))/(max(abs(GN)));
            %     plot(impz(num{test},denom{test},length(GN)));
            
        end
        
%                 figure
%                 hold on
%                 plot(GN)
%                 plot(GI)

        % make sure there is still relevant signal
        if max(abs(GN))<.15
            adequateRMSE = 1;
        end
        %if there is at least one RMSE that is below .15, then continue,
        %otherwise, remove one sample from the beginning of the signal.
        if any(RMSE<.15) 
          adequateRMSE = 1;
        else
            GN = GN(2:end);
            removed = removed +1;
                            
        end
    end
    
    %% 4) find the best impulse model
    %compare the RMSE to find the best model with a minimum improvement of
    %5%
    best = 1;
    for p = 2:4
        change = RMSE(best)-RMSE(p);
        if change>.05
            best =  p;
        end
    end
    [~,Poles] = prony(GN,0,best);
    Poles = Poles(2:end);
    RMax = max(Poles);
    %make sure that the pole is less than .89
    if RMax>.89
        isPSO = 0;
        REnd = 'highPole';
        continue
    end
    % Check that the amplitude is large enough to be a PSO
    Amp = max(abs(GI(:,best)));
    if Amp<.15
        isPSO = 0;
        REnd = 'lowAmp';
        continue
    end
   %% 5) Find the point where the signal decays below threshold
   % this was originally based on a calculation that used the decay of the
   % signal, but now simply has a threshold which oscillations with an
   % amplitude below the threshold are considered to likely be noise.
   
    signal = GI(:,best);
    %gets all the points that are outside the threshold of .08
    belowThresh = [0;abs(signal)<ThreshAngle;0]; % add 0's for future searches
    changes = find(belowThresh==0);
    %finds the point where the signal was below the threshold for at least
    %6 ms, that is the endpoint
   
    if~isempty(changes)
        DipIndicies = [find(diff(changes)>3);length(changes)];
        DipsBelowSignal = changes(DipIndicies(1));
        
    else  DipsBelowSignal = length(signal); 
        % if no period below threshold for long enough, then leave it empty
    end
       
    % if the signal is always above threshold, then make the end point the
    % last point of the signal
    endpoint = min(DipsBelowSignal,length(signal));
    
    % add in the  points that were previously removed
    endpoint = removed + endpoint ;
    % if there has been no endpoint found, set it as the length of the
    % sample
    if isempty(endpoint);
        endpoint = length(GN);
    end
    if length(Time)<endpoint
        endpoint = length(Time);
        pause
    end
    %% 6) Calculate the amplitude duration ratio
    %Now calculate the ratio of the amplitudes, and make sure that it is
    %not just a slow movement, by adding up the max positive and negative
    %amplitudes, and dividing by the time.
    S =  (max(PPSO(1:endpoint)) - min(PPSO(1:endpoint))) ./ (endpoint*.002);

    if S< 15
        isPSO = 2;
        continue
    end
    PSOendpoint = Time(endpoint);
    isPSO = 2;
    
    %     %if S is less than 400/2*Sampling Frequency, (in this case it is probably a slow
    %     %movement
    %     200/500
    %     figure
    %     hold on
    %     plot(GN)
    %     plot(decay)
    %     plot(GI(:,best))
    
end

end


