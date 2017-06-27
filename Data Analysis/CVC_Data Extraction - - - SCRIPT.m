

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
% - This first line is to select the folder I want
folder_name = uigetdir

% - Now that you've selected the folder, execute all these lines of code
A = folder_name;
pth = strcat(A,'\');
rawlist = getrawdata(pth,'randycode');
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


Var = rawlist.data;

[VALUE_DATA] = Con_Value_Choice_Data(Var);

%Var1 = c_rawlist;
TS = VALUE_DATA.trial_structure;
[CVC_TEMPORAL_DATA] = CVC_Temportal_Dat(Var,TS);
[CVC_BOUT_PAUSE_DATA] = CVC_Bout_Pause_Dat(Var,TS);


