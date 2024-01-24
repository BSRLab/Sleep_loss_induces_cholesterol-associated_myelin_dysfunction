%% SLEEP ANALYSIS
clear all

%%folder
outpath= 'C:\Users\Utente\Desktop';
inpath= 'C:\Users\Utente\Desktop';

ID= 'Animal';
STRAIN= 'W';
MANIPULATION= 'BASELINE';

%choose time bin 
bin= 24;
epoch= 4;
fs= 512;
consolidate= 2;
offset_onset= 0;

% upload file hypnogram from scoring
%anal_pathin= [inpath filesep STRAIN filesep MANIPULATION filesep ID filesep 'MAT_files'];
anal_pathin=inpath;
file_hypno= fullfile(anal_pathin,['hypnogram_',ID,'_',MANIPULATION,'.mat']);
load(file_hypno);

all_lenght=24*3600/4;
blocks=bin*3600/epoch;
nblocks=length(hypno)/blocks;
sz=[blocks nblocks];
states=unique(hypno);
nstates=length(states);
%reshape hypno
hypno_blocked=reshape(hypno,sz);
%number of epoch scored with each state
h=1:nblocks;
for j=1:nblocks
    blockvar(j)=strcat("Block_",num2str(j));
end
epoch_state=zeros(nstates,nblocks);
for i=1:nstates
            epoch_state(i, :) = sum(hypno_blocked==states(i));
end

Vigilance_State_all= string(states);
epoch_state=table(Vigilance_State_all,epoch_state);
epoch_state = splitvars(epoch_state,'epoch_state','NewVariableNames',blockvar);

%table for seconds in each state
seconds_state=epoch_state;
seconds_state{:,2:end}= seconds_state{:,2:end}*epoch;

%table for percentage in each state
perc_state=epoch_state;
perc_state{:,2:end}= perc_state{:,2:end}*100./repmat(sum(epoch_state{:,2:end}),(length(states)),1);

%merging states with artifacts in order to have REM, NREM and WAKE
W=epoch_state.Vigilance_State_all =="W"| epoch_state.Vigilance_State_all=="WA";
NR=epoch_state.Vigilance_State_all =="NR"| epoch_state.Vigilance_State_all=="NA";
R=epoch_state.Vigilance_State_all =="R"| epoch_state.Vigilance_State_all=="RA";
UN=epoch_state.Vigilance_State_all =="UN";
statewoa={'NREM','REM','WAKE','UNSCORED'};
Vigilance_State=string(statewoa)';
wake=sum(epoch_state{W,2:end},1);
nrem=sum(epoch_state{NR,2:end},1);
rem=sum(epoch_state{R,2:end},1);
un= sum(epoch_state{UN,2:end},1);
%epochs in each state
epoch_state_all=table(Vigilance_State,[nrem;rem;wake;un]);
epoch_state_all= splitvars(epoch_state_all,'Var2','NewVariableNames',blockvar);

seconds_state_all=epoch_state_all;
seconds_state_all{:,2:end}= seconds_state_all{:,2:end}*epoch;

perc_state_all=epoch_state_all;
perc_state_all{:,2:end}=perc_state_all{:,2:end}*100./repmat(wake+nrem+rem+un,length(statewoa),1);


hypno_number=zeros(length(hypno),1);
for z=1:length(hypno)
if hypno(z)=='WA'|hypno(z)=='W';
    hypno_number(z)=3;
elseif hypno(z)=='NA'|hypno(z)=='NR';
    hypno_number(z)=1;
elseif hypno(z)=='RA'|hypno(z)=='R';
    hypno_number(z)=2;
elseif hypno(z)=='UN';
    hypno_number(z)=4;
else
    hypno_number(z)=NaN;
end
end

%remove NaN
hypno_number(isnan(hypno_number))=[];

%Vector with consecutive occurrence counts of the states
consec_occurence= diff([0;find(diff(hypno_number));numel(hypno_number)]);
consec_state=zeros(length(consec_occurence),3);
position=0;
for j=1:length(consec_occurence)
    consec_state(j,1)=consec_occurence(j);
    position=position+consec_occurence(j);
    consec_state(j,2)=hypno_number(position);
    consec_state(j,3)=position-consec_occurence(j)+1;
end


for h=1:length(statewoa)
state_num_episodes_all(h,1)=length(find(consec_state(:,2)==h&consec_state(:,1)>=consolidate));
state_episodes_mean_duration_all(h,1)= mean(consec_state(consec_state(:,2)==h&consec_state(:,1)>=consolidate,1));
state_episodes_std_duration_all(h,1)= std(consec_state(consec_state(:,2)==h&consec_state(:,1)>=consolidate,1));
end

state_num_episodes_all=table(Vigilance_State,state_num_episodes_all);

state_episodes_mean_duration_all=table(Vigilance_State,state_episodes_mean_duration_all);

state_episodes_std_duration_all=table(Vigilance_State,state_episodes_std_duration_all);

%NUMBER OF EPISODES AND DURATION in TIME BIN
rv = cell(length(statewoa), sz(2)); 
for l=1:length(consec_state)
    c_stage=consec_state(l,2);
    if (consec_state(l,3)/blocks)==sz(2)
        c_block= floor((consec_state(l,3)/blocks));
    else
    c_block = floor((consec_state(l,3)/blocks))+1;
    end 
    c_len=consec_state(l,1);
    rv{c_stage, c_block} = [rv{c_stage, c_block}; c_len];        
end
    
state_num_episodes=zeros(length(statewoa), sz(2));
state_mean_episodes_duration=zeros(length(statewoa), sz(2));
state_std_epidodes_duration=zeros(length(statewoa), sz(2));

for a=1:length(statewoa)
    for b=1:sz(2)
        logic=rv{a,b}>=consolidate;
        state_num_episodes(a,b)=sum(logic);
        state_mean_episodes_duration(a,b)=mean(rv{a,b}(logic)*epoch);%in seocnds
        state_std_epidodes_duration(a,b)=std(rv{a,b}(logic)*epoch);
    end
end

%table for the number of episodes
state_num_episodes=table(Vigilance_State,state_num_episodes);
state_num_episodes= splitvars(state_num_episodes,'state_num_episodes','NewVariableNames',blockvar);

%table for mean durations
state_mean_episodes_duration=table(Vigilance_State,state_mean_episodes_duration);
state_mean_episodes_duration= splitvars(state_mean_episodes_duration,'state_mean_episodes_duration','NewVariableNames',blockvar);

%table for standard deviation of duration
state_std_epidodes_duration=table(Vigilance_State,state_std_epidodes_duration);
state_std_epidodes_duration= splitvars(state_std_epidodes_duration,'state_std_epidodes_duration','NewVariableNames',blockvar);



%%OUTPUT 
%epochs for each state
SUMMARY.QUANTITY.epoch_state= epoch_state;
%seconds for each state
SUMMARY.QUANTITY.seconds_state= seconds_state; 
%percentage for each state
SUMMARY.QUANTITY.perc_state= perc_state;
%epochs for each state (merging artifacts)
SUMMARY.QUANTITY.epoch_state_all= epoch_state_all;
%seconds for each state (merging artifacts)
SUMMARY.QUANTITY.seconds_state_all=seconds_state_all;
%percentage for each state (merging artifacts)
SUMMARY.QUANTITY.perc_state_all= perc_state_all;
% number of episodes
SUMMARY.ARCHITECTURE.state_num_episodes= state_num_episodes;
%mean duration of episodes
SUMMARY.ARCHITECTURE.state_mean_episodes_duration= state_mean_episodes_duration; 
%std of episodes duration
SUMMARY.ARCHITECTURE.state_std_episodes_duration=state_std_epidodes_duration;
%number of episodes for each total state
SUMMARY.ARCHITECTURE.state_num_episodes_all=state_num_episodes_all;
%mean duration for each total state 
SUMMARY.ARCHITECTURE.state_episodes_mean_duration_all=state_episodes_mean_duration_all;
%std of episodes duration for each total state
SUMMARY.ARCHITECTURE.state_episodes_std_duration_all=state_episodes_std_duration_all;

%save file
outpath_mat=[outpath filesep  ['BIN_', num2str(bin), '_hr'] filesep MANIPULATION filesep STRAIN filesep ID filesep 'MAT_FILES'];
outpath_xls=[outpath filesep  ['BIN_', num2str(bin), '_hr'] filesep MANIPULATION filesep STRAIN filesep ID filesep 'EXCEL_FILES'];
file_mat=fullfile(outpath_mat,['SUMMARY_','BIN_', num2str(bin), '_hr_',ID,'_',STRAIN,'_',MANIPULATION,'.mat']);
if not(isfolder(outpath_mat))
    mkdir(outpath_mat)
end
if not(isfolder(outpath_xls))
    mkdir(outpath_xls)
end

%save .mat file
save(file_mat,'SUMMARY')

%QUANTITY EXCEL
file_quantity_xls=fullfile(outpath_xls,['SLEEP_QUANTITY_','BIN_', num2str(bin), '_hr_',ID,'_',STRAIN,'_',MANIPULATION,'.xlsx']);
q_f=fields(SUMMARY.QUANTITY);

for xx=1:length(q_f)
    warning('off','MATLAB:xlswrite:AddSheet')
    writetable(SUMMARY.QUANTITY.(q_f{xx}),file_quantity_xls,'Sheet',q_f{xx})
end

%ARCHITECTURE EXCEL
file_arch_xls=fullfile(outpath_xls,['SLEEP_ARCHITECTURE_','BIN_', num2str(bin), '_hr_',ID,'_',STRAIN,'_',MANIPULATION,'.xlsx']);
a_f=fields(SUMMARY.ARCHITECTURE);
xl_n=a_f;
for yy=1:length(a_f)-1
    warning('off','MATLAB:xlswrite:AddSheet')
    g_lent=length(a_f{yy})-31;
    if g_lent>0
        xl_n{yy}(end-g_lent)=[];
    end
    writetable(SUMMARY.ARCHITECTURE.(a_f{yy}),file_arch_xls,'Sheet',xl_n{yy})
end

