function [DATA] = CVC_Bout_Pause_Dat(Var,TS)
%UNTITLED2 Summary of this function goes here

% - Var = rawlist.data
% - TS = trial_structure from VALUE_DATA


%   Detailed explanation goes here
% - - This is a function which will determine the duration of every lever
% press, the inter-response-times (IRT's) or Pauses, as well as the
% structure of the bouts of responding: Number of Bouts, Lenth of each Bout
% Ect


% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
n = length(Var);

%LLeverOn = 0027 \***
%RLeverOn= 0028  \***

%ForcedTrial = 3300
%ChoiceTrial = 3301

%PelletTrial = 3320
%SucroseTrial = 3321

%Pelletchoice = 0100\code for hold choice
%Sucrosechoice = 0101\code for press choice
%ForcedOpOut = 03330\opting out of forced trial

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


% - - - All just shit to allow me to know dimensions of TS - - - 
%=======================================================%
%tells how many start trial/end trials are there%
for i = 1:n
    num_start{1,i} = length(find(Var{1,i}(:,2)==111));
    num_end{1,i} = length(find(Var{1,i}(:,2)==112));
end

%Defines rows in Ru of Start, end, AND End Session%
for i = 1:n
    start_rows{1,i} = find(Var{1,i}(:,2)==111);
    end_rows{1,i} = find(Var{1,i}(:,2)==112);
    end_all{1,i} = find(Var{1,i}(:,2)==118);
end

%Abbrevating rows with start trial%
sr = start_rows;


%Subtracting 1 from start trial b/c sr-1 is end of trial%
for i = 1:n
    sr_m1{1,i} = sr{1,i}(:,1) - 1;
end

%removing first row%
for i = 1:n
    sr_m1{1,i}(1,:)=[];
end


%Adding position of session end to end of sr_m1, b/c thats the
%row # for the end of the last trial%
for i = 1:n
    sr_m1{1,i}(end+1,:) = end_all{1,i};
end

%Variable defining # of start and end trials%
for i = 1:n
    Sess_Info_num_trials{1,i} = [num_start{1,i};num_end{1,i}];
end
SI_nt = Sess_Info_num_trials;   %easier var name to work with%

%Easy name for how many trials are there%
for i = 1:n
    A{1,i} = SI_nt{1,i}(1,1);
end

%Making a TRIALS data structure: each row contains all time stamps
%and event codes from that trial number%

%NOTE: switched "i" here and now "j" = # subs (i.e.13)
for j = 1:n
    for i = 1:A{1,j}
        trial_structure{i,j} = Var{1,j}(sr{1,j}(i,1):sr_m1{1,j}(i,1),:);
    end
end
% - - - All just shit to allow me to know dimensions of TS - - -

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

TS_Lever_Zero = TS;

%Easier name to use
TSL = TS_Lever_Zero;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -
%Seting the start time of each trial to zero

            % Two Steps: A. and B. 

% A. 
%Figuring out the start time of each trial (Globally)
for i = 1:n
    for j = 1:A{1,i}
        TSL_Trial_Start_Global_Time(j,i) =  TSL{j,i}(1,1);
    end
end 

% B.
%Subtracting this time from all of the other times in the TSL struct within
%each trial - Zeroing everything so each trial begins at time zero

for i = 1:n
    for j = 1:A{1,i}
        TSL{j,i}(:,1) = TSL{j,i}(:,1) - TSL_Trial_Start_Global_Time(j,i);
    end
end 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -


% - Determining the row number that the lever comes out for each trial
for i = 1:n
    for j = 1:A{1,i}
        if length(find(TSL{j,i}(:,2) == 27)) == 0
            TSL_LeverOut_Rows{j,i} = find(TSL{j,i}(:,2) == 28);
        else
            TSL_LeverOut_Rows{j,i} = find(TSL{j,i}(:,2) == 27);
        end
    end
end

% - Determining the time that the lever comes out for each trial

for i = 1:n
    for j = 1:A{1,i}
        TSL_LeverOut_Time(j,i) = TSL{j,i}(TSL_LeverOut_Rows{j,i},1);
    end
end


% - -  Getting just the forced trials
for i = 1:n
    for j = 1:10
        TSL_Forced{j,i} = TSL{j,i};
    end
end

TSLF = TSL_Forced;

% - Determining the Press Rows
for i = 1:n
    for j = 1:10
        if length(find(TSL{j,i}(:,2) == 28)) == 0
            TSL_Any_Press_Rows{j,i} = find(TSLF{j,i}(:,2) == 1015);
        else
            TSL_Any_Press_Rows{j,i} = find(TSLF{j,i}(:,2) == 1016);
        end
    end
end



for i = 1:n
    P = 1;
    H = 1;
    for j = 1:10
        if length(find(TSLF{j,i}(:,2) == 3320)) == 0
            SUCROSE_STRUCT{P,i} = TSLF{j,i};
            P = P + 1;
        else
            PELLET_STRUCT{H,i} = TSLF{j,i};
            H = H + 1;
        end
    end
end

% - - Determining the Pellet Rows and Hold Rows

for i = 1:n
    for j = 1:5
        if length(find(PELLET_STRUCT{j,i}(:,2) == 28)) == 0
            Pellet_Down_Rows{j,i} = find(PELLET_STRUCT{j,i}(:,2) == 1015);
            Pellet_Up_Rows{j,i} = find(PELLET_STRUCT{j,i}(:,2) == 1017);
            Forced_Pellet_Lever_Out_Rows(j,i) = find(PELLET_STRUCT{j,i}(:,2) == 27);
        else
            Pellet_Down_Rows{j,i} = find(PELLET_STRUCT{j,i}(:,2) == 1016);
            Pellet_Up_Rows{j,i} = find(PELLET_STRUCT{j,i}(:,2) == 1018);
            Forced_Pellet_Lever_Out_Rows(j,i) = find(PELLET_STRUCT{j,i}(:,2) == 28);
        end
    end
end


for i = 1:n
    for j = 1:5
        if length(find(SUCROSE_STRUCT{j,i}(:,2) == 28)) == 0
            Sucrose_Down_Rows{j,i} = find(SUCROSE_STRUCT{j,i}(:,2) == 1015);
            Sucrose_Up_Rows{j,i} = find(SUCROSE_STRUCT{j,i}(:,2) == 1017);
            Forced_Sucrose_Lever_Out_Rows(j,i) = find(SUCROSE_STRUCT{j,i}(:,2) == 27);
        else
            Sucrose_Down_Rows{j,i} = find(SUCROSE_STRUCT{j,i}(:,2) == 1016);
            Sucrose_Up_Rows{j,i} = find(SUCROSE_STRUCT{j,i}(:,2) == 1018);
            Forced_Sucrose_Lever_Out_Rows(j,i) = find(SUCROSE_STRUCT{j,i}(:,2) == 28);
        end
    end
end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

for i = 1:n
    for j = 1:5
        PELLET_DOWN_TIMES{j,i} = PELLET_STRUCT{j,i}(Pellet_Down_Rows{j,i},1);
        PELLET_UP_TIMES{j,i} = PELLET_STRUCT{j,i}(Pellet_Up_Rows{j,i},1);
        SUCROSE_Down_TIMES{j,i} = SUCROSE_STRUCT{j,i}(Sucrose_Down_Rows{j,i},1);
        SUCROSE_Up_TIMES{j,i} = SUCROSE_STRUCT{j,i}(Sucrose_Up_Rows{j,i},1);
        Pellet_num(j,i) = length(Pellet_Down_Rows{j,i});
        Sucrose_num(j,i) = length(Sucrose_Down_Rows{j,i});
        Pellet_Lever_Out_Time(j,i) = PELLET_STRUCT{j,i}(Forced_Pellet_Lever_Out_Rows(j,i),1);
        Sucrose_Lever_Out_Time(j,i) = SUCROSE_STRUCT{j,i}(Forced_Sucrose_Lever_Out_Rows(j,i),1);
    end
end

for i = 1:n
    for j = 1:5
        Sucrose_num1(j,i) = length(Sucrose_Up_Rows{j,i});
        Pellet_num(j,i) = length(Pellet_Up_Rows{j,i});
    end
end

for i = 1:n
    for j = 1:5
       if length(find(PELLET_STRUCT{j,i}(:,2) == 3330) == 1)
            PELLET_DURATIONS{j,i} = NaN;
       else    
            PELLET_DURATIONS{j,i} = PELLET_UP_TIMES{j,i} - PELLET_DOWN_TIMES{j,i};
       end
    end
end

for i = 1:n
    PELLET_DURATION_DATA{1,i} = [PELLET_DURATIONS{1,i};PELLET_DURATIONS{2,i};PELLET_DURATIONS{3,i};...
        PELLET_DURATIONS{4,i};PELLET_DURATIONS{5,i}];
end


for i = 1:n
   Mean_Pellet_Response_Durations(1,i) = nanmean(PELLET_DURATION_DATA{1,i}); 
end


for i = 1:n
   x_dur_pell{1,i} = 1:1:length(PELLET_DURATION_DATA{1,i}); 
end

for i = 1:n
    for j = 1:5
       if length(find(SUCROSE_STRUCT{j,i}(:,2) == 3330) == 1)
            SUCROSE_DURATIONS{j,i} = NaN;
       else    
            SUCROSE_DURATIONS{j,i} = SUCROSE_Up_TIMES{j,i} - SUCROSE_Down_TIMES{j,i};
       end
    end
end

for i = 1:n
    SUCROSE_DURATION_DATA{1,i} = [SUCROSE_DURATIONS{1,i};SUCROSE_DURATIONS{2,i};SUCROSE_DURATIONS{3,i};...
        SUCROSE_DURATIONS{4,i};SUCROSE_DURATIONS{5,i}];
end


for i = 1:n
   Mean_Sucrose_Response_Durations(1,i) = nanmean(SUCROSE_DURATION_DATA{1,i}); 
end


% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
%- - Need Actual Press Up/Down Times to get IRT's
% - Three steps:
    % A) Getting the actual Times
    % B) Removing the first Down Time and Last UP Time
    % C) Difference between UP(n) and Down(n - 1)

    
% = = = For Pellets = = =     
for i = 1:n
    for j = 1:5
       if length(find(PELLET_STRUCT{j,i}(:,2) == 3330) == 1)
            PELLET_FOR_IRT_UP_TIMES{j,i} = NaN;
            PELLET_FOR_IRT_DOWN_TIMES{j,i} = NaN;
       else    
            PELLET_FOR_IRT_UP_TIMES{j,i} = PELLET_UP_TIMES{j,i};
            PELLET_FOR_IRT_DOWN_TIMES{j,i} = PELLET_DOWN_TIMES{j,i};
       end
    end
end

for i = 1:n
    for j = 1:5
       if length(find(PELLET_STRUCT{j,i}(:,2) == 3330) == 1)
            PELLET_FOR_IRT_UP_TIMES{j,i} = NaN;
            PELLET_FOR_IRT_DOWN_TIMES{j,i} = NaN;
       else    
            PELLET_FOR_IRT_UP_TIMES{j,i}(end,:) = [];
            PELLET_FOR_IRT_DOWN_TIMES{j,i}(1,:) = [];
       end
    end
end


for i = 1:n
    for j = 1:5
       if length(find(PELLET_STRUCT{j,i}(:,2) == 3330) == 1)
            PELLET_IRTs{j,i} = NaN;
       else    
            PELLET_IRTs{j,i} = PELLET_FOR_IRT_DOWN_TIMES{j,i} - PELLET_FOR_IRT_UP_TIMES{j,i};
       end
    end
end

for i = 1:n
    PELLET_IRT_DATA{1,i} = [PELLET_IRTs{1,i};PELLET_IRTs{2,i};PELLET_IRTs{3,i};...
        PELLET_IRTs{4,i};PELLET_IRTs{5,i}];
end

for i = 1:n
   x_irt_pell{1,i} = 1:1:length(PELLET_IRT_DATA{1,i}); 
end

aa = 1:1:9;
aa = aa';
x_irt_num = [aa;aa;aa;aa;aa];



% = = = For Sucrose = = =     
for i = 1:n
    for j = 1:5
       if length(find(SUCROSE_STRUCT{j,i}(:,2) == 3330) == 1)
            SUCROSE_FOR_IRT_UP_TIMES{j,i} = NaN;
            SUCROSE_FOR_IRT_DOWN_TIMES{j,i} = NaN;
       else    
            SUCROSE_FOR_IRT_UP_TIMES{j,i} = SUCROSE_Up_TIMES{j,i};
            SUCROSE_FOR_IRT_DOWN_TIMES{j,i} = SUCROSE_Down_TIMES{j,i};
       end
    end
end

for i = 1:n
    for j = 1:5
       if length(find(SUCROSE_STRUCT{j,i}(:,2) == 3330) == 1)
            SUCROSE_FOR_IRT_UP_TIMES{j,i} = NaN;
            SUCROSE_FOR_IRT_DOWN_TIMES{j,i} = NaN;
       else    
            SUCROSE_FOR_IRT_UP_TIMES{j,i}(end,:) = [];
            SUCROSE_FOR_IRT_DOWN_TIMES{j,i}(1,:) = [];
       end
    end
end


for i = 1:n
    for j = 1:5
       if length(find(SUCROSE_STRUCT{j,i}(:,2) == 3330) == 1)
            SUCROSE_IRTs{j,i} = NaN;
       else    
            SUCROSE_IRTs{j,i} = SUCROSE_FOR_IRT_DOWN_TIMES{j,i} - SUCROSE_FOR_IRT_UP_TIMES{j,i};
       end
    end
end

for i = 1:n
    SUCROSE_IRT_DATA{1,i} = [SUCROSE_IRTs{1,i};SUCROSE_IRTs{2,i};SUCROSE_IRTs{3,i};...
        SUCROSE_IRTs{4,i};SUCROSE_IRTs{5,i}];
end

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

% - - Finding Pauses and Bout Lengths


% = = = = NUMBER OF PAUSES AND PAUSE DURATIONS = = = = = = = 

% - Structure of Pause Rows for Different Durations
% - Row Storage Schedule: 1sec, 2sec, 5sec, 10sec
for i = 1:n
    Pellet_Pause_rows{1,i} = find(PELLET_IRT_DATA{1,i}(:,1) >= 1);
    Sucrose_Pause_rows{1,i} = find(SUCROSE_IRT_DATA{1,i}(:,1) >= 1);
    
    Pellet_Pause_rows{2,i} = find(PELLET_IRT_DATA{1,i}(:,1) >= 2);
    Sucrose_Pause_rows{2,i} = find(SUCROSE_IRT_DATA{1,i}(:,1) >= 2);
    
    Pellet_Pause_rows{3,i} = find(PELLET_IRT_DATA{1,i}(:,1) >= 5);
    Sucrose_Pause_rows{3,i} = find(SUCROSE_IRT_DATA{1,i}(:,1) >= 5);
    
    Pellet_Pause_rows{4,i} = find(PELLET_IRT_DATA{1,i}(:,1) >= 10);
    Sucrose_Pause_rows{4,i} = find(SUCROSE_IRT_DATA{1,i}(:,1) >= 10);
end

% =========================================================================================
% =========================================================================================
% = = Getting IRT's after removing the pauses 

Filtered_Pellet_IRT_DATA = {PELLET_IRT_DATA;PELLET_IRT_DATA;...
    PELLET_IRT_DATA;PELLET_IRT_DATA};
Filtered_Sucrose_IRT_DATA = {SUCROSE_IRT_DATA;SUCROSE_IRT_DATA;...
    SUCROSE_IRT_DATA;SUCROSE_IRT_DATA};

% = Getting IRT Data and Actually Removing the Different Pause Def's
for i = 1:n
    for j = 1:4
        if isempty(Pellet_Pause_rows{j,i})
            Filtered_Pellet_IRT_DATA{j,1}{1,i} = Filtered_Pellet_IRT_DATA{j,1}{1,i};
        else
            Filtered_Pellet_IRT_DATA{j,1}{1,i}(Pellet_Pause_rows{j,i}(:,1),:) = [];
        end
    end
end

for i = 1:n
    for j = 1:4
        if isempty(Sucrose_Pause_rows{j,i})
            Filtered_Sucrose_IRT_DATA{j,1}{1,i} = Filtered_Sucrose_IRT_DATA{j,1}{1,i};
        else
            Filtered_Sucrose_IRT_DATA{j,1}{1,i}(Sucrose_Pause_rows{j,i}(:,1),:) = [];
        end
    end
end


for i = 1:n
    for j = 1:4
        Mean_Pellet_IRT(j,i) = nanmean(Filtered_Pellet_IRT_DATA{j,1}{1,i});
    end
end

for i = 1:n
    for j = 1:4
        Mean_Sucrose_IRT(j,i) = nanmean(Filtered_Sucrose_IRT_DATA{j,1}{1,i});
    end
end
% =========================================================================================
% =========================================================================================



% - Determining the number of Pauses made
% - Row Storage Schedule: 1sec, 2sec, 5sec, 10sec
for i = 1:n
    for j = 1:4
        Pellet_Pause_num{j,i} = length(Pellet_Pause_rows{j,i});
        Sucrose_Pause_num{j,i} = length(Sucrose_Pause_rows{j,i});
    end
end

% - - Creating a Variable of just the "Pause" Durations
for i = 1:n
    for j = 1:4
        Pellet_Pause_Durations{j,i} = PELLET_IRT_DATA{1,i}(Pellet_Pause_rows{j,i},1);
        Sucrose_Pause_Durations{j,i} = SUCROSE_IRT_DATA{1,i}(Sucrose_Pause_rows{j,i},1);
    end
end

for i = 1:n
    for j = 1:4
        Mean_Pellet_Pause_Duration{j,i} = nanmean(Pellet_Pause_Durations{j,i});
        Mean_Sucrose_Pause_Duration{j,i} = nanmean(Sucrose_Pause_Durations{j,i});
        Sum_Pellet_Pause_Duration{j,i} = sum(Pellet_Pause_Durations{j,i});
        Sum_Sucrose_Pause_Duration{j,i} = sum(Sucrose_Pause_Durations{j,i});
    end
end

% = = = = = = Bout Lengths = = = = = = = = = 
% - Row Storage Schedule: 1sec, 2sec, 5sec, 10sec
for i = 1:n
    for j = 1:5

            Trial_Pellet_Pause_rows{j,i}{1,1} = find(PELLET_IRTs{j,i}(:,1) >= 1);
            Trial_Sucrose_Pause_rows{j,i}{1,1} = find(SUCROSE_IRTs{j,i}(:,1) >= 1);
            
            Trial_Pellet_Pause_rows{j,i}{2,1} = find(PELLET_IRTs{j,i}(:,1) >= 2);
            Trial_Sucrose_Pause_rows{j,i}{2,1} = find(SUCROSE_IRTs{j,i}(:,1) >= 2);
            
            Trial_Pellet_Pause_rows{j,i}{3,1} = find(PELLET_IRTs{j,i}(:,1) >= 5);
            Trial_Sucrose_Pause_rows{j,i}{3,1} = find(SUCROSE_IRTs{j,i}(:,1) >= 5);
            
            Trial_Pellet_Pause_rows{j,i}{4,1} = find(PELLET_IRTs{j,i}(:,1) >= 10);
            Trial_Sucrose_Pause_rows{j,i}{4,1} = find(SUCROSE_IRTs{j,i}(:,1) >= 10);
            
            Pel_No_Pause_Detector{j,i} = find(PELLET_IRTs{j,i}(:,1) >= .00001);
    end
end


% = = = TWO STEPS
% - A) Getting Opt Out Trials to Contain NaN's
% - B) Getting Trials with no Pauses to be able to lead to proper
% calculation of the bout lengths

% - A)
% - Getting trials which were optouts to contain NaN's
for i = 1:n
    for j = 1:5
        for k = 1:4
            if isnan(PELLET_DURATIONS{j,i})
                Trial_Pellet_Pause_rows{j,i}{k,1} = NaN;
            else
                Trial_Pellet_Pause_rows{j,i}{k,1} = Trial_Pellet_Pause_rows{j,i}{k,1};
            end
        end
    end
end

%- B) Getting Trials with no pauses to be able to calculate correct bout
%lengths

% - Number of Pauses per Trial
for i = 1:n
    for j = 1:5
        for k = 1:4    
            if isnan(Trial_Pellet_Pause_rows{j,i}{k,1})
                trial_num_pellet_pause{j,i}{k,1} = NaN;
            else
                trial_num_pellet_pause{j,i}{k,1} = length(Trial_Pellet_Pause_rows{j,i}{k,1});
            end
        end
    end
end


for i = 1:n
    for j = 1:5
        for k = 1:4    
            if isnan(SUCROSE_DURATIONS{j,i})
                trial_num_sucrose_pause{j,i}{k,1} = NaN;
            else
                trial_num_sucrose_pause{j,i}{k,1} = length(Trial_Sucrose_Pause_rows{j,i}{k,1});
            end
        end
    end
end

for i = 1:n
    for j = 1:5
        for k = 1:4    
            if isnan(Trial_Sucrose_Pause_rows{j,i}{k,1})
                trial_num_sucrose_pause{j,i}{k,1} = NaN;
            else
                trial_num_sucrose_pause{j,i}{k,1} = length(Trial_Sucrose_Pause_rows{j,i}{k,1});
            end
        end
    end
end

% - - Preparing the data so I can extract the bout lengths from it
% - - 2 Step Process
        % A) Adding 1 to the pause rows so I can just subtract them and get the
               % bout lengths
        % B) Adding the total number of responses to the end 
        
% - - Getting the number needed for part B
for i = 1:n
    for j = 1:5
        for k = 1:4    
            if isnan(PELLET_DURATIONS{j,i})
                Num_Pell_Presses{j,i} = NaN;
            else
                Num_Pell_Presses{j,i} = length(PELLET_DURATIONS{j,i});
            end
         end
     end
end
 
for i = 1:n
    for j = 1:5
        for k = 1:4    
            if isnan(SUCROSE_DURATIONS{j,i})
                Num_Suc_Presses{j,i} = NaN;
            else
                Num_Suc_Presses{j,i} = length(SUCROSE_DURATIONS{j,i});
            end
         end
     end
 end

% - Implementing Part A and B
for i = 1:n
    for j = 1:5
        for k = 1:4    
            if isnan(Num_Pell_Presses{j,i})
                Pellet_bout_length_trials_PRE{j,i}{k,1} = NaN;
            else
                Pellet_bout_length_trials_PRE{j,i}{k,1} = [1;Trial_Pellet_Pause_rows{j,i}{k,1};Num_Pell_Presses{j,i}];
            end
        end
    end
end

for i = 1:n
    for j = 1:5
        for k = 1:4    
            if isnan(Trial_Sucrose_Pause_rows{j,i}{k,1})
                Sucrose_bout_length_trials_PRE{j,i}{k,1} = NaN;
            else
                Sucrose_bout_length_trials_PRE{j,i}{k,1} = [1;Trial_Sucrose_Pause_rows{j,i}{k,1};Num_Suc_Presses{j,i}];
            end
        end
    end
end


% - - - Determining the length of the bouts within each trial
for i = 1:n
    for j = 1:5
        for k = 1:4    
            if isnan(Trial_Pellet_Pause_rows{j,i}{k,1})
                Pellet_Bout_Length_Trials{j,i}{k,1} = NaN;
            else
                Pellet_Bout_Length_Trials{j,i}{k,1} = diff(Pellet_bout_length_trials_PRE{j,i}{k,1});
            end
        end
    end
end

for i = 1:n
    for j = 1:5
        for k = 1:4    
            if isnan(Trial_Sucrose_Pause_rows{j,i}{k,1})
                Sucrose_Bout_Length_Trials{j,i}{k,1} = NaN;
            else
                Sucrose_Bout_Length_Trials{j,i}{k,1} = diff(Sucrose_bout_length_trials_PRE{j,i}{k,1});
            end
        end
    end
end

% - Making variable with all pellet bout legnths over all trials
for i = 1:n
    for k = 1:4
        PELLET_BOUT_LENGTH_DATA{k,i} = [Pellet_Bout_Length_Trials{1,i}{k,1};...
            Pellet_Bout_Length_Trials{2,i}{k,1};...
            Pellet_Bout_Length_Trials{3,i}{k,1};...
            Pellet_Bout_Length_Trials{4,i}{k,1};...
            Pellet_Bout_Length_Trials{5,i}{k,1}];
    end
end



% - - Calculating the average pellet bout length
for i = 1:n
        Var_BOUT_1sec{1,i} = PELLET_BOUT_LENGTH_DATA{1,i};
        Var_BOUT_2sec{1,i} = PELLET_BOUT_LENGTH_DATA{2,i};
        Var_BOUT_5sec{1,i} = PELLET_BOUT_LENGTH_DATA{3,i};
        Var_BOUT_10sec{1,i} = PELLET_BOUT_LENGTH_DATA{4,i};
end

[Pel_BOUT_DF_1sec] = trials_cellmat(Var_BOUT_1sec,Var);
[Pel_BOUT_DF_2sec] = trials_cellmat(Var_BOUT_2sec,Var);
[Pel_BOUT_DF_5sec] = trials_cellmat(Var_BOUT_5sec,Var);
[Pel_BOUT_DF_10sec] = trials_cellmat(Var_BOUT_10sec,Var);


pb_DF_1sec = Pel_BOUT_DF_1sec.Trials_Matrix;
pb_DF_2sec = Pel_BOUT_DF_2sec.Trials_Matrix;
pb_DF_5sec = Pel_BOUT_DF_5sec.Trials_Matrix;
pb_DF_10sec = Pel_BOUT_DF_10sec.Trials_Matrix;


for i = 1:n
    Pel_NUM_BOUT(1,i) = length(pb_DF_1sec(:,i));
    Pel_NUM_BOUT(2,i) = length(pb_DF_2sec(:,i));
    Pel_NUM_BOUT(3,i) = length(pb_DF_5sec(:,i));
    Pel_NUM_BOUT(4,i) = length(pb_DF_10sec(:,i));
end

for i = 1:n
    Bout_NaN{1,i} = length(find(isnan(pb_DF_1sec(:,i)) == 1));
    Bout_NaN{2,i} = length(find(isnan(pb_DF_2sec(:,i)) == 1));
    Bout_NaN{3,i} = length(find(isnan(pb_DF_5sec(:,i)) == 1));
    Bout_NaN{4,i} = length(find(isnan(pb_DF_10sec(:,i)) == 1));
end

for i = 1:n
    Pellet_Number_of_BOUTS(1,i) = Pel_NUM_BOUT(1,i) - Bout_NaN{1,i};
    Pellet_Number_of_BOUTS(2,i) = Pel_NUM_BOUT(2,i) - Bout_NaN{2,i};
    Pellet_Number_of_BOUTS(3,i) = Pel_NUM_BOUT(3,i) - Bout_NaN{3,i};
    Pellet_Number_of_BOUTS(4,i) = Pel_NUM_BOUT(4,i) - Bout_NaN{4,i};
end

Mean_Pellet_Bout_Length(1,:) = nanmean(pb_DF_1sec);
Mean_Pellet_Bout_Length(2,:) = nanmean(pb_DF_2sec);
Mean_Pellet_Bout_Length(3,:) = nanmean(pb_DF_5sec);
Mean_Pellet_Bout_Length(4,:) = nanmean(pb_DF_10sec);


% - Making a vriable with all sucrose bout lengths combined
% - Making a vriable with all sucrose bout lengths combined
for i = 1:n
    for k = 1:4
        SUCROSE_BOUT_LENGTH_DATA{k,i} = [Sucrose_Bout_Length_Trials{1,i}{k,1};...
            Sucrose_Bout_Length_Trials{2,i}{k,1};...
            Sucrose_Bout_Length_Trials{3,i}{k,1};...
            Sucrose_Bout_Length_Trials{4,i}{k,1};...
            Sucrose_Bout_Length_Trials{5,i}{k,1}];
    end
end

% - Determining the average sucrose bout length
for i = 1:n
    Var_BOUT_suc_1sec{1,i} = SUCROSE_BOUT_LENGTH_DATA{1,i};
    Var_BOUT_suc_2sec{1,i} = SUCROSE_BOUT_LENGTH_DATA{2,i};
    Var_BOUT_suc_5sec{1,i} = SUCROSE_BOUT_LENGTH_DATA{3,i};
    Var_BOUT_suc_10sec{1,i} = SUCROSE_BOUT_LENGTH_DATA{4,i};
end

[Suc_BOUT_DF_1sec] = trials_cellmat(Var_BOUT_suc_1sec,Var);
[Suc_BOUT_DF_2sec] = trials_cellmat(Var_BOUT_suc_2sec,Var);
[Suc_BOUT_DF_5sec] = trials_cellmat(Var_BOUT_suc_5sec,Var);
[Suc_BOUT_DF_10sec] = trials_cellmat(Var_BOUT_suc_10sec,Var);

sb_DF_1sec = Suc_BOUT_DF_1sec.Trials_Matrix;
sb_DF_2sec = Suc_BOUT_DF_2sec.Trials_Matrix;
sb_DF_5sec = Suc_BOUT_DF_5sec.Trials_Matrix;
sb_DF_10sec = Suc_BOUT_DF_10sec.Trials_Matrix;

for i = 1:n
    Suc_NUM_BOUT(1,i) = length(sb_DF_1sec(:,i));
    Suc_NUM_BOUT(2,i) = length(sb_DF_2sec(:,i));
    Suc_NUM_BOUT(3,i) = length(sb_DF_5sec(:,i));
    Suc_NUM_BOUT(4,i) = length(sb_DF_10sec(:,i));
end

for i = 1:n
    Bout_NaN_suc{1,i} = length(find(isnan(sb_DF_1sec(:,i)) == 1));
    Bout_NaN_suc{2,i} = length(find(isnan(sb_DF_2sec(:,i)) == 1));
    Bout_NaN_suc{3,i} = length(find(isnan(sb_DF_5sec(:,i)) == 1));
    Bout_NaN_suc{4,i} = length(find(isnan(sb_DF_10sec(:,i)) == 1));
end

for i = 1:n
    Sucrose_Number_of_BOUTS(1,i) = Suc_NUM_BOUT(1,i) - Bout_NaN_suc{1,i};
    Sucrose_Number_of_BOUTS(2,i) = Suc_NUM_BOUT(2,i) - Bout_NaN_suc{2,i};
    Sucrose_Number_of_BOUTS(3,i) = Suc_NUM_BOUT(3,i) - Bout_NaN_suc{3,i};
    Sucrose_Number_of_BOUTS(4,i) = Suc_NUM_BOUT(4,i) - Bout_NaN_suc{4,i};
end

Mean_Sucrose_Bout_Length(1,:) = nanmean(sb_DF_1sec);
Mean_Sucrose_Bout_Length(2,:) = nanmean(sb_DF_2sec);
Mean_Sucrose_Bout_Length(3,:) = nanmean(sb_DF_5sec);
Mean_Sucrose_Bout_Length(4,:) = nanmean(sb_DF_10sec);




% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
for i = 1:n
    for j = 1:5
        if length(find(PELLET_STRUCT{j,i}(:,2) == 22)) == 1
            PELLET_Reward_Rows(j,i) = find(PELLET_STRUCT{j,i}(:,2) == 22);
        
        else
            PELLET_Reward_Rows(j,i) = NaN;
        end
    end
end

for i = 1:n
    for j = 1:5
        if length(find(SUCROSE_STRUCT{j,i}(:,2) == 25)) == 1
            SUCROSE_Reward_Rows(j,i) = find(SUCROSE_STRUCT{j,i}(:,2) == 25);
        else
            SUCROSE_Reward_Rows(j,i) = NaN;
        end
    end
end


for i = 1:n
    for j = 1:5
        if length(find(PELLET_STRUCT{j,i}(:,2) == 22)) == 1
            PELLET_Reward_Times(j,i) = PELLET_STRUCT{j,i}(PELLET_Reward_Rows(j,i),1); 
        else
            PELLET_Reward_Times(j,i) = NaN;
        end
    end
end

for i = 1:n
    for j = 1:5
        if length(find(SUCROSE_STRUCT{j,i}(:,2) == 25)) == 1
            SUCROSE_Reward_Times(j,i) = SUCROSE_STRUCT{j,i}(SUCROSE_Reward_Rows(j,i),1);
        else
            SUCROSE_Reward_Times(j,i) = NaN;
        end
    end
end



% - - Determining Times to reward from:
% - - 1) First press to reward
% - - 2) Trial start to reward

for i = 1:n
   for j = 1:5
      if isempty(PELLET_DOWN_TIMES{j,i}) == 1 
          Working_Pellet_Time(j,i) =  NaN;
          Total_Pellet_Time(j,i) =  NaN;
      else
          Working_Pellet_Time(j,i) =  PELLET_Reward_Times(j,i) - PELLET_DOWN_TIMES{j,i}(1,1);
          Total_Pellet_Time(j,i) =  PELLET_Reward_Times(j,i) - Pellet_Lever_Out_Time(j,i);
      end
   end
end

for i = 1:n
   for j = 1:5
       if isempty(SUCROSE_Down_TIMES{j,i}) == 1
           Working_Sucrose_Time(j,i) =  NaN;
           Total_Sucrose_Time(j,i) =  NaN;
       else
            Working_Sucrose_Time(j,i) =  SUCROSE_Reward_Times(j,i) - SUCROSE_Down_TIMES{j,i}(1,1);
            Total_Sucrose_Time(j,i) =  SUCROSE_Reward_Times(j,i) - Sucrose_Lever_Out_Time(j,i);
       end
   end
end


% - - Determining Latency to start working

for i = 1:n
   for j = 1:5
        if isempty(PELLET_DOWN_TIMES{j,i}) == 1 
              Latency_to_First_Pellet(j,i) = NaN;
        else
              Latency_to_First_Pellet(j,i) = PELLET_DOWN_TIMES{j,i}(1,1) - Pellet_Lever_Out_Time(j,i);
        end
   end
end

for i = 1:n
   for j = 1:5
      if isempty(SUCROSE_Down_TIMES{j,i}) == 1
            Latency_to_First_Sucrose(j,i) =  NaN;
      else    
            Latency_to_First_Sucrose(j,i) =  SUCROSE_Down_TIMES{j,i}(1,1) - Sucrose_Lever_Out_Time(j,i);
      end
   end
end


for i = 1:n
   Sub_PELLET_Reward_Times{1,i} = PELLET_Reward_Times(:,i); 
   Sub_SUCROSE_Reward_Times{1,i} = SUCROSE_Reward_Times(:,i);
end

SPRT = Sub_PELLET_Reward_Times;
SHRT = Sub_SUCROSE_Reward_Times;

pn = 5;
a = 1:1:pn;
a = a';

for i = 1:n
    Full_Pellet_Trials{1,i} = find(Pellet_num(:,i) == 20);
end

fpt = Full_Pellet_Trials;

for i = 1:n
    lfp(1,i) = length(Full_Pellet_Trials{1,i});
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - Getting Averages over Data

for i = 1:16
    Mean_Working_Pellet_Time = nanmean(Working_Pellet_Time);
    Mean_Total_Pellet_Time = nanmean(Total_Pellet_Time);
    Mean_Latency_to_First_Pellet = nanmean(Latency_to_First_Pellet);
    
    Mean_Working_Sucrose_Time = nanmean(Working_Sucrose_Time);
    Mean_Total_Sucrose_Time = nanmean(Total_Sucrose_Time);
    Mean_Latency_to_First_Sucrose = nanmean(Latency_to_First_Sucrose);
end

% - - - - - - - - - - DATA OUTPUT OF THE FUNCTION - - - - - - - - - - - %


DATA = struct('PELLET_DURATIONS_DATA',{{}},...
    'SUCROSE_DURATIONS_DATA',{{}},...
    'Mean_Pell_DURATIONS',{{}},...
    'Mean_Suc_DURATIONS',{{}},...
    'Pell_IRTs',{{}},...
    'Suc_IRTs',{{}},...
    'Mean_Pell_IRTs',{{}},...
    'Mean_Suc_IRTs',{{}},...
    'Pellet_Pause_Durations',{{}},...
    'Sucrose_Pause_Durations',{{}},...
    'Pellet_Pause_num',{{}},...
    'Sucrose_Pause_num',{{}},...
    'Mean_Pellet_Pause_Duration',{{}},...
    'Mean_Sucrose_Pause_Duration',{{}},...
    'Sum_Pellet_Pause_Duration',{{}},...
    'Sum_Sucrose_Pause_Duration',{{}},...
    'trial_num_pellet_pause',{{}},...
    'trial_num_sucrose_pause',{{}},...
    'Pellet_Bout_Length_Trials',{{}},...
    'Sucrose_Bout_Length_Trials',{{}},...
    'PELLET_BOUT_LENGTH_DATA',{{}},...
    'SUCROSE_BOUT_LENGTH_DATA',{{}},...
    'Mean_Pellet_Bout_Length',{{}},...
    'Mean_Sucrose_Bout_Length',{{}}); 


DATA.PELLET_DURATIONS_DATA = PELLET_DURATION_DATA;   %Pellet Response Duations
DATA.SUCROSE_DURATIONS_DATA	= SUCROSE_DURATION_DATA; %Sucrose Response Durations

DATA.Mean_Pell_DURATIONS = Mean_Pellet_Response_Durations;   %Mean response duration
DATA.Mean_Suc_DURATIONS	= Mean_Sucrose_Response_Durations;

DATA.Pell_IRTs = Filtered_Pellet_IRT_DATA; %List of All Pellet IRT's
DATA.Suc_IRTs = Filtered_Sucrose_IRT_DATA;
DATA.Mean_Pell_IRTs = Mean_Pellet_IRT; %Mean of All IRT's
DATA.Mean_Suc_IRTs = Mean_Sucrose_IRT;

DATA.Pellet_Pause_Durations	= Pellet_Pause_Durations; %Trial x trial - Pause Durations
DATA.Sucrose_Pause_Durations = Sucrose_Pause_Durations;

DATA.Pellet_Pause_num = Pellet_Pause_num;
DATA.Sucrose_Pause_num = Sucrose_Pause_num;

DATA.Mean_Pellet_Pause_Duration = Mean_Pellet_Pause_Duration; %Trial x Trial
DATA.Mean_Sucrose_Pause_Duration = Mean_Sucrose_Pause_Duration;	%Trial x Trial
DATA.Sum_Pellet_Pause_Duration = Sum_Pellet_Pause_Duration;	%Trial x Trial
DATA.Sum_Sucrose_Pause_Duration = Sum_Sucrose_Pause_Duration;	%Trial x Trial


DATA.trial_num_pellet_pause	= trial_num_pellet_pause; %Trial x Trial # Pauses
DATA.trial_num_sucrose_pause = trial_num_sucrose_pause;	%Trial x Trial # Pauses

DATA.Pellet_Bout_Length_Trials = Pellet_Bout_Length_Trials;	%Trial x Trial bout length		
DATA.Sucrose_Bout_Length_Trials = Sucrose_Bout_Length_Trials;	%Trial x Trial bout length


DATA.PELLET_BOUT_LENGTH_DATA = PELLET_BOUT_LENGTH_DATA;	%All bout lengths combined
DATA.SUCROSE_BOUT_LENGTH_DATA = SUCROSE_BOUT_LENGTH_DATA;	%%All bout lengths combined


DATA.Mean_Pellet_Bout_Length =	Mean_Pellet_Bout_Length;	% Averaged over all trials
DATA.Mean_Sucrose_Bout_Length = Mean_Sucrose_Bout_Length;

end

