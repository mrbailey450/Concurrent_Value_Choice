function [Con_Value_DATA] = Con_Value_Choice_Data(Var)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


n = length(Var);

% - - - Overall Session Variables - - - %

% Number of Reinforcers%
for i = 1:n
    Sucrose_Reins(1,i) = length(find(Var{1,i}(:,2) == 25));
end;

for i = 1:n
    Pellet_Reins(1,i) = length(find(Var{1,i}(:,2) == 22));
end;

for i = 1:n
   Total_Reinforcers(1,i) = Sucrose_Reins(1,i) + Pellet_Reins(1,i); 
end

%Number of OpOuts
for i = 1:n
    OpOuts(1,i) = length(find(Var{1,i}(:,2) == 3330));
end;

%Session Duration%
for i = 1:n;
    start(1,i) = find(Var{1,i}(:,2) == 113);
    stop(1,i) = find(Var{1,i}(:,2) == 114);
    SessionDur(1,i) = (Var{1,i}(stop(1,i),1) - Var{1,i}(start(1,i),1))/60;
end;


% ====================================================================%
% =====================Getting the Trial Structure====================%



%Determining the rows of each trial ending time, bc it's 
%always preceded by a lever up from a succ press
for j = 1:n
    end_trial{1,j} = find(Var{1,j}(:,2)==112);
end

e = end_trial;

%Getting variable with Lengths of end trial variable
%to use in for command as the larger value

for j = 1:n
    l(1,j) = length(e{1,j}(:,1));
end

%Making the "e minus 1" variable

for j = 1:n
for i = 1:l(1,j)
    if e{1,j}(i,1) >= 0.1
        em{1,j}(i,1) = e{1,j}(i,1) - 1;
    end
end
end


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

% ====================================================================%
% ====================================================================%
% ====================================================================%

Var1 = trial_structure;


for i = 1:n 
    num_trials(1,i) = length(find(Var{1,i}(:,2) == 111));
end

l = num_trials;

for i = 1:n
    for j = 1:l(1,i)
        if isempty(Var1{j,i})
            Var1{j,i} = [1,1];
        else
             Var1{j,i} = Var1{j,i};
        end
    end        
end    

for i = 1:n
    for j = 1:l(1,i)
        Pellet_row{j,i} = find(Var1{j,i}(:,2)==100);
        Sucrose_row{j,i} = find(Var1{j,i}(:,2)==25);
    end
end


for i = 1:n
    for j = 1:l(1,i)
        if isempty(Pellet_row{j,i})
            Pellet_trial(j,i) = 0;
        else
            Pellet_trial(j,i) = 1;
        end
    end        
end    


for i = 1:n
    for j = 1:l(1,i)
        if isempty(Sucrose_row{j,i})
            Sucrose_trial(j,i) = 0;
        else
            Sucrose_trial(j,i) = 1;
        end
    end        
end  


for i = 1:n
    Session_Pellet(1,i) = length(find(Pellet_trial(:,i)==1));
    Session_Sucrose(1,i) = length(find(Sucrose_trial(:,i)==1));
end
        
for i = 1:n
    Cum_Sess_Pellet(:,i) = cumsum(Pellet_trial(:,i));
    Cum_Sess_Sucrose(:,i) = cumsum(Sucrose_trial(:,i));
end        

% %Pellets in last 20 Trials - - Choice Trials%
% for i = 1:n
%     if Total_Reinforcers(1,i) >= 30
%         Pellets_in_Choice(1,i) = sum(Pellet_trial(end-20:end,i));
%     else
%         Pellets_in_Choice(1,i) = sum(Pellet_trial(1 + OpOuts(1,i):end,i));
%     end
% end  


for i = 1:n
    for j = 1:l(1,i)
        Choice_row{j,i} = find(Var1{j,i}(:,2)==3301);
    end
end


for i = 1:n
    for j = 1:l(1,i)
        if isempty(Choice_row{j,i})
            Choice_trial(j,i) = 0;
        else    
            Choice_trial(j,i) = 1;
        end
    end        
end

ll = max(length(Choice_trial));

for i = 1:n
    for j = 1:ll
        if isempty(trial_structure{j,i})
            Choice_trial(j,i) = NaN;
        else    
            Choice_trial(j,i) = Choice_trial(j,i);
        end
    end        
end


for i = 1:n
    TEST{1,i} = find(Choice_trial(:,i) == 1);
end

for i = 1:n
    Size_Choice(:,i) = size(TEST{1,i});
end

for i = 1:n
   Total_Choice(1,i) = Size_Choice(1,i); 
end


for i = 1:n
       Pellets_in_Choice(1,i) = sum(Pellet_trial(TEST{1,i},i));
       
end


for i = 1:n
    Sucrose_in_Choice(1,i) = sum(Sucrose_trial(TEST{1,i},i));
end

for i = 1:n
    Prop_Sucrose_Num(1,i) = Sucrose_in_Choice(1,i) ./ Total_Choice(1,i);
end

%Start row for each trial%     ----NOTE: A{1,j}-1 bc we don't care about the
                                    %last trial
for j = 1:n
    for i = 1: A{1,j}-1;
        start_row{1,j}(i,1) = find(Var1{i,j}(:,2)==111);
    end
end

for j = 1:n
    for i = 1: A{1,j}-1;
        end_row{1,j}(i,1) = find(Var1{i,j}(:,2)==121);
    end
end


%Time to complete each trial% 
for j = 1:n
    for i = 1: A{1,j}-1;
        Time_to_Complete{i,j} = trial_structure{i,j}(end_row{1,j}(i,1),1)-...
            trial_structure{i,j}(start_row{1,j}(i,1),1);
    end
end


lt = max(length(Time_to_Complete));

for j = 1:n
    for i = 1:lt;
        if isempty(Time_to_Complete{i,j})
            Time_to_Complete{i,j} = NaN;
        else
            Time_to_Complete{i,j} = Time_to_Complete{i,j};
        end
    end        
end 

TIME_TO_COMPLETE = cell2mat(Time_to_Complete);

%Var1 = trial_structure;
%x = rawlist.data;
%n = length(x);

l = Total_Reinforcers;

for i = 1:n
    for j = 1:l
        if isempty(Var1{j,i})
            Var1{j,i} = [1,1];
        else
             Var1{j,i} = Var1{j,i};
        end
    end        
end    

for i = 1:n
    for j = 1:l
        Pellet_row{j,i} = find(Var1{j,i}(:,2)==100);
    end
end



for i = 1:n
    Session_Pellet(1,i) = length(find(Pellet_trial(:,i)==1));
end
        
for i = 1:n
    Cum_Sess_Pellet(:,i) = cumsum(Pellet_trial(:,i));
end        
       


for i = 1:n
    for j = 1:l
        L_press_rows{j,i} = find(Var1{j,i}(:,2)==1015);
        R_press_rows{j,i} = find(Var1{j,i}(:,2)==1016);
    end
end

for i = 1:n
    for j = 1:l
        rein_row{j,i} = find(Var1{j,i}(:,2)==25);
    end
end


for i = 1:n
    for j = 1:l
        Left_Sizes{j,i} = size(L_press_rows{j,i}); 
   end
end

for i = 1:n
    for j = 1:l
        Right_Sizes{j,i} = size(R_press_rows{j,i}); 
   end
end


for i = 1:n
    for j = 1:l
        if Left_Sizes{j,i}(1,1)>Right_Sizes{j,i}(1,1)
            Active_Press_Rows{j,i} = L_press_rows{j,i};
        else
            Active_Press_Rows{j,i} = R_press_rows{j,i};
        end
    end
end
APR = Active_Press_Rows;

for i = 1:n
    for j = 1:l
        no_trial(j,i) = isempty(rein_row{j,i});
    end
end


for i = 1:n
    for j = 1:l
        if no_trial(j,i) ==1;
            time_between{j,i} = NaN;
        else
            time_between{j,i} = Var1{j,i}(rein_row{j,i},1) - Var1{j,i}(APR{j,i}(1,1),1); 
        end
    end
end

TIME_BETWEEN = cell2mat(time_between);

for i = 1:n
   CUm_TIme_Between(:,i) = cumsum(TIME_BETWEEN(:,i));
end

% Determining number of reins and opt outs%%

for i = 1:n
    for j = 1:l
        Rein_row{j,i} = find(Var1{j,i}(:,2)==25);
    end
end


for i = 1:n
    for j = 1:l
        if isempty(Rein_row{j,i})
            Rein_trial(j,i) = 0;
        else
             Rein_trial(j,i) = 1;
        end
    end        
end   


%Reins in last 30%
for i = 1:n
       Rein_30(1,i) = sum(Rein_trial(11:end,i));
end

for i = 1:n
    Choice(1,i) = length(find(Var{1,i}(:,2) == 3301));
end

for i = 1:n
   Pellet_Choice_NUM(1,i) = sum(Pellet_trial(11:30,i));
end



% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% - - - Output of the Function - - - - 

Con_Value_DATA = struct('Prop_Sucrose_Num',[],'OpOuts',[],...
    'Total_Reinforcers',[],'Sucrose_Reins',[],'Pellet_Reins',[],...
    'SessionDur',[],'TIME_TO_COMPLETE',[],'trial_structure',{{}});

Con_Value_DATA.Prop_Sucrose_Num = Prop_Sucrose_Num;
Con_Value_DATA.OpOuts = OpOuts;
Con_Value_DATA.Total_Reinforcers = Total_Reinforcers;
Con_Value_DATA.Sucrose_Reins = Sucrose_Reins;
Con_Value_DATA.Pellet_Reins = Pellet_Reins;
Con_Value_DATA.SessionDur = SessionDur;
Con_Value_DATA.TIME_TO_COMPLETE = TIME_TO_COMPLETE;
Con_Value_DATA.trial_structure = trial_structure;

%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

end

