function [CVC_TEMPORAL_DATA] = CVC_Temportal_Dat(Var,TS)
%UNTITLED2 Summary of this function goes here

%   Var = rawlist.data
%   TS = VALUE_DATA.trial_structure
%   

%   Detailed explanation goes here

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

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
    for j = 1:5
       if length(find(SUCROSE_STRUCT{j,i}(:,2) == 3330) == 1)
            SUCROSE_DURATIONS{j,i} = NaN;
       else    
            SUCROSE_DURATIONS{j,i} = SUCROSE_Up_TIMES{j,i} - SUCROSE_Down_TIMES{j,i};
       end
    end
end


% Pellet_down_for_IRT = PELLET_DOWN_TIMES;
% Pellet_up_for_IRT = PELLET_UP_TIMES;
% 
% pd_IRT = Pellet_down_for_IRT;
% pu_IRT = Pellet_up_for_IRT;
% 
% for i = 1:n
%     for j = 1:5
%         if isempty(pd_IRT{j,i})
%             pd_IRT{j,i} = NaN;
%         else
%             pd_IRT{j,i}(1,:) = [];
%         end
%     end
% end
% 
% for i = 1:n
%     for j = 1:5
%         if isempty(pu_IRT{j,i})
%             pu_IRT{j,i} = NaN;
%         else
%             pu_IRT{j,i}(end,:) = [];
%         end    
%     end
% end
% 
% 
% for i = 1:n
%     for j = 1:5
%         if isempty(pu_IRT{j,i})
%             IRT_Pellet_Presses{j,i} = NaN;
%         else
%             IRT_Pellet_Presses{j,i} = pd_IRT{j,i} - pu_IRT{j,i};
%         end
%     end
% end



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



% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
% - - - - - - - - - - - Output of the Function - - - - - - - - - - - - 

CVC_TEMPORAL_DATA = struct('PELLET_DOWN_TIMES',{{}},...
    'SUCROSE_Down_TIMES',{{}},...
    'SUCROSE_Up_TIMES',{{}},...
    'Pellet_num',[],...
    'Sucrose_num',[],...
    'SUCROSE_DURATIONS',{{}},...
    'PELLET_Reward_Times',[],...
    'SUCROSE_Reward_Times',[],...
    'Working_Pellet_Time',[],...
    'Total_Pellet_Time',[],...
    'Working_Sucrose_Time',[],...
    'Total_Sucrose_Time',[],...
    'Latency_to_First_Pellet',[],...
    'Latency_to_First_Sucrose',[],...
    'Data_Structs',{{}},...
    'Means',{{}});


CVC_TEMPORAL_DATA.PELLET_DOWN_TIMES = PELLET_DOWN_TIMES;
CVC_TEMPORAL_DATA.SUCROSE_Down_TIMES = SUCROSE_Down_TIMES;
CVC_TEMPORAL_DATA.SUCROSE_Up_TIMES = SUCROSE_Up_TIMES;
CVC_TEMPORAL_DATA.Pellet_num = Pellet_num;
CVC_TEMPORAL_DATA.Sucrose_num = Sucrose_num;

CVC_TEMPORAL_DATA.SUCROSE_DURATIONS = SUCROSE_DURATIONS;


CVC_TEMPORAL_DATA.PELLET_Reward_Times = PELLET_Reward_Times;
CVC_TEMPORAL_DATA.SUCROSE_Reward_Times = SUCROSE_Reward_Times;


CVC_TEMPORAL_DATA.Working_Pellet_Time = Working_Pellet_Time;
CVC_TEMPORAL_DATA.Total_Pellet_Time = Total_Pellet_Time;
CVC_TEMPORAL_DATA.Working_Sucrose_Time = Working_Sucrose_Time;
CVC_TEMPORAL_DATA.Total_Sucrose_Time = Total_Sucrose_Time;


CVC_TEMPORAL_DATA.Latency_to_First_Pellet = Latency_to_First_Pellet;
CVC_TEMPORAL_DATA.Latency_to_First_Sucrose = Latency_to_First_Sucrose;

% - - The Pellet and Sucrose Data Structs for Forced Trials - - - - 
CVC_TEMPORAL_DATA.Data_Structs = struct('PELLET_STRUCT',{{}},'SUCROSE_STRUCT',{{}});
CVC_TEMPORAL_DATA.Data_Structs.PELLET_STRUCT = PELLET_STRUCT;
CVC_TEMPORAL_DATA.Data_Structs.SUCROSE_STRUCT = SUCROSE_STRUCT;


CVC_TEMPORAL_DATA.Means = struct('Mean_Working_Pellet_Time',[],...
    'Mean_Total_Pellet_Time',[],...
    'Mean_Latency_to_First_Pellet',[],...
    'Mean_Working_Sucrose_Time',[],...
    'Mean_Total_Sucrose_Time',[],...
    'Mean_Latency_to_First_Sucrose',[]);




CVC_TEMPORAL_DATA.Means.Mean_Working_Pellet_Time = Mean_Working_Pellet_Time;
CVC_TEMPORAL_DATA.Means.Mean_Total_Pellet_Time = Mean_Total_Pellet_Time;
CVC_TEMPORAL_DATA.Means.Mean_Latency_to_First_Pellet = Mean_Latency_to_First_Pellet;
    
CVC_TEMPORAL_DATA.Means.Mean_Working_Sucrose_Time = Mean_Working_Sucrose_Time;
CVC_TEMPORAL_DATA.Means.Mean_Total_Sucrose_Time = Mean_Total_Sucrose_Time;
CVC_TEMPORAL_DATA.Means.Mean_Latency_to_First_Sucrose = Mean_Latency_to_First_Sucrose;

end

