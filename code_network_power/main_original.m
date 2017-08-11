clc;clear all;
%cvx_solver sedumi
%cvx_solver sdpt3
%cvx_solver mosek
cvx_solver sdpt3
%cvx_solver_settings('max_iters', 10^5, 'scale', 1,'eps',1.00e-5);
%cvx_solver_settings('SCALE', 1);
cvx_quiet(true)

%% Problem Data (I)
channel_generating=0;
LC=10; %# loops for channel realizarions

CBF=1;
IRLS=1;
IRLS_2=1;
IRLS_3=0;
GSBF_Alternating=0;
Baseline=1;
Exhaustive=1;


%%
Area=4*10^3;
L=12; %Number of RAUs
N1=2;  % Antennas of each RAU

M=5; %Number of Multicast Groups
K1=2; %Number of Mobile Users in Each Multicast Group

%Q=[6:-2:0]';  %QoS in dB
Q=[8:-2:0]';

%%
N_set=N1*ones(L,1); % Set of Antennas for all the RAU
N=sum(N_set);
K_set=K1*ones(M,1); %Set of Numbers of Mobile Users in the Multicast Groups
K=sum(K_set);


Pc=3*ones(L,1)+[0:1:(L-1)]';  %power consumption of the fronthaul link
%Pc=5.6*ones(L,1);
eita=1/4; %amplifier efficiency coefficient


%% Problem Data for Params (II)
params.L=L;
params.N=N;
params.M=M;
params.K=K;
params.K_set=K_set;
params.N_set=N_set;
params.delta_set=ones(M,1);
params.P_set=10^(0)*ones(L,1);

%%%%%%initialize results%%%%%%%%%%%%%%%%%%%%%%
TotalPower_CBF_temp=0; TotalPower_GSBFAlternating_temp=0; TotalPower_IRLSAlternating_temp=0;
TotalPower_Baseline_temp=0; TotalPower_Exhaustive_temp=0; %network power consumption
TotalPower_IRLSAlternating_2_temp=0;TotalPower_IRLSAlternating_3_temp=0;

TransmitPower_CBF_temp=0;  TransmitPower_GSBFAlternating_temp=0;TransmitPower_IRLSAlternating_temp=0;
 TransmitPower_Baseline_temp=0; TransmitPower_Exhaustive_temp=0; %Total Transmit Power consumption
TransmitPower_IRLSAlternating_2_temp=0;TransmitPower_IRLSAlternating_3_temp=0;
 
A_number_Heuristic_temp=0;  A_number_GSBFAlternating_temp=0; A_number_IRLSAlternating_temp=0;
A_number_Baseline_temp=0; A_number_Exhaustive_temp=0;  %Optimal number of active RAUs
A_number_IRLSAlternating_2_temp=0;A_number_IRLSAlternating_3_temp=0;
               
               
SU_counter=0;
%% params for IRLS
p1=1;p2=0.5;p3=0.05;e=0.001;


%% Channel Generation
if channel_generating==true
    %% fix large scale
for ss=1:LC   %generate channel
    for m=1:length(K_set)
        D(:,m,:)=largefading(Area,L,K_set(m),N_set);
    end
    for m=1:1:length(K_set)
        H(:,m,:,ss)=smallfading(Area,L,K_set(m),N_set, D(:,m,:));
    end
end
%save('H.mat','H');
%% standard channel
% for lc=1:1:LC
%     for m=1:1:length(K_set)
%         H(:,m,:,lc)=channel_realization(Area,L,K_set(m),N_set);
%     end
% end
% save('H.mat','H');
end
load('H.mat');
H=H(:,:,:,1:10);

for lp=1:LC
    
    SU_channel_temp=1;  %recode if the channel is feasible
    
    params.H=H(:,:,:,lp);  %NxMxK channel matrix
    
    for lq=1:length(Q)
        params.r_set=10^(Q(lq)/10)*ones(M,1);  %set of SINR thresholds  for each group 
%% Coordinated Beamforming%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CBF==true

%% Problem Solving
params.Inactive_index=[]; 
params.Active_index=[1:L];
params.weight=(1/eita)*ones(length(params.Active_index),1);
params.rankone=true; %return rankone solution
[V_CBF, feasible_CBF] = powermin_cvx(params);

%% Output for the Particular QoS
if feasible_CBF==1
TotalPower_CBF(lq)=(1/eita)*norm(V_CBF,'fro')^2+sum(Pc);
TransmitPower_CBF(lq)=(1/eita)*norm(V_CBF,'fro')^2;
else
    SU_channel_temp=0;  % this channel if infeasible
    break;
end


end 

%% IRLS algorithm wiht p=1
if IRLS==1

alternating_max=20; alternating_abs=10^99;
alternating=0;
value1=10^99; value2=0;
w=1*ones(L,1);
    while alternating<=alternating_max
        
        alternating=alternating+1; %interation numbers
        value1=value2; %recode the old objective value
        
        params.Inactive_index=[]; 
        params.Active_index=[1:L];
        params.weight=w;
        params.rankone=false; %return general rank solution
        [Q_IRLS, feasible_IRLS] = powermin_cvx(params);
        
        if feasible_IRLS==false
            break;
        else
            value_temp=zeros(L,1);
            for l=1:L 
                for m=1:M
                    value_temp(l)=value_temp(l)+trace(Q_IRLS(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
                end                
            end
            %value2=params.w'*value_temp;
            %value2=sum(value_temp);
            value2=sum((value_temp+e^2).^(p1/2));
            w=Pc.*(value_temp+e^2).^(p1/2-1);
            alternating_abs(alternating)=abs(value1-value2); %absolute value of the adjacent objective values
        end
    end
    
    if feasible_IRLS==false
        SU_channel_temp=0;  %This channel is infeasible
        break;
    else  
        %% RRH Ordering Based on the Approximated Group Sparse Beamformer
        Value_IRLS=zeros(L,1); 
        for l=1:L 
            temp1=0;
            for m=1:M
                temp1=temp1+trace(Q_IRLS(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
            end
            Value_IRLS(l)=(Pc(l))^(-1)*temp1;
        end
        
        %% Sparsity Parttern Based RRH Ordering 
      [BS,BS_index]=sort(Value_IRLS);
 
        %%
        D_set=[]; A_set=BS_index;  %A_set: active RRH set, D_set: inactive RRH set
        
        %% Process Deflation Procedure
        for A=0:L-1
            
            params.Inactive_index=D_set; 
            params.Active_index=A_set;
            params.weight=(1/eita)*ones(L,1);
            params.rankone=true; %return rankone solution
            [Wsolution, feasible] = powermin_cvx(params); %power minimization given the active RRH set
            
            if feasible==1
                
                TotalPower_IRLSAlternating(lq)=(1/eita)*norm(Wsolution,'fro')^2+sum(Pc(A_set));  %recode current values
                TransmitPower_IRLSAlternating(lq)=(1/eita)*norm(Wsolution,'fro')^2;
                A_number_IRLSAlternating(lq)=length(A_set);
                
                D_set=[D_set, BS_index(A+1)];  %%%process next selection: updata active RRH set and number
                A_set=setdiff([1:L],D_set); 
                
            else
                break;
            end
        end   
    end  
end

%% IRLS algorithm wiht p=0.5
if IRLS_2==1

alternating_max=20; alternating_abs=10^99;
alternating=0;
value1=10^99; value2=0;
w=1*ones(L,1);
    while alternating<=alternating_max
        
        alternating=alternating+1; %interation numbers
        value1=value2; %recode the old objective value
        
        params.Inactive_index=[]; 
        params.Active_index=[1:L];
        params.weight=w;
        params.rankone=false; %return general rank solution
        [Q_IRLS, feasible_IRLS] = powermin_cvx(params);
        
        if feasible_IRLS==false
            break;
        else
            value_temp=zeros(L,1);
            for l=1:L 
                for m=1:M
                    value_temp(l)=value_temp(l)+trace(Q_IRLS(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
                end                
            end
            %value2=params.w'*value_temp;
            %value2=sum(value_temp);
            value2=sum((value_temp+e^2).^(p2/2));
            w=Pc.*(value_temp+e^2).^(p2/2-1);
            alternating_abs(alternating)=abs(value1-value2); %absolute value of the adjacent objective values
        end
    end
    
    if feasible_IRLS==false
        SU_channel_temp=0;  %This channel is infeasible
        break;
    else  
        %% RRH Ordering Based on the Approximated Group Sparse Beamformer
        Value_IRLS=zeros(L,1); 
        for l=1:L 
            temp1=0;
            for m=1:M
                temp1=temp1+trace(Q_IRLS(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
            end
            Value_IRLS(l)=(Pc(l))^(-1)*temp1;
        end
        
        %% Sparsity Parttern Based RRH Ordering 
      [BS,BS_index]=sort(Value_IRLS);
 
        %%
        D_set=[]; A_set=BS_index;  %A_set: active RRH set, D_set: inactive RRH set
        
        %% Process Deflation Procedure
        for A=0:L-1
            
            params.Inactive_index=D_set; 
            params.Active_index=A_set;
            params.weight=(1/eita)*ones(L,1);
            params.rankone=true; %return rankone solution
            [Wsolution, feasible] = powermin_cvx(params); %power minimization given the active RRH set
            
            if feasible==1
                
                TotalPower_IRLSAlternating_2(lq)=(1/eita)*norm(Wsolution,'fro')^2+sum(Pc(A_set));  %recode current values
                TransmitPower_IRLSAlternating_2(lq)=(1/eita)*norm(Wsolution,'fro')^2;
                A_number_IRLSAlternating_2(lq)=length(A_set);
                
                D_set=[D_set, BS_index(A+1)];  %%%process next selection: updata active RRH set and number
                A_set=setdiff([1:L],D_set); 
                
            else
                break;
            end
        end   
    end  
end 

%% IRLS algorithm wiht p=0.05
if IRLS_3==1

alternating_max=20; alternating_abs=10^99;
alternating=0;
value1=10^99; value2=0;
w=1*ones(L,1);
    while alternating<=alternating_max
        
        alternating=alternating+1; %interation numbers
        value1=value2; %recode the old objective value
        
        params.Inactive_index=[]; 
        params.Active_index=[1:L];
        params.weight=w;
        params.rankone=false; %return general rank solution
        [Q_IRLS, feasible_IRLS] = powermin_cvx(params);
        
        if feasible_IRLS==false
            break;
        else
            value_temp=zeros(L,1);
            for l=1:L 
                for m=1:M
                    value_temp(l)=value_temp(l)+trace(Q_IRLS(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
                end                
            end
            %value2=params.w'*value_temp;
            %value2=sum(value_temp);
            value2=sum((value_temp+e^2).^(p3/2));
            w=Pc.*(value_temp+e^2).^(p3/2-1);
            alternating_abs(alternating)=abs(value1-value2); %absolute value of the adjacent objective values
        end
    end
    
    if feasible_IRLS==false
        SU_channel_temp=0;  %This channel is infeasible
        break;
    else  
        %% RRH Ordering Based on the Approximated Group Sparse Beamformer
        Value_IRLS=zeros(L,1); 
        for l=1:L 
            temp1=0;
            for m=1:M
                temp1=temp1+trace(Q_IRLS(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
            end
            Value_IRLS(l)=(Pc(l))^(-1)*temp1;
        end
        
        %% Sparsity Parttern Based RRH Ordering 
      [BS,BS_index]=sort(Value_IRLS);
 
        %%
        D_set=[]; A_set=BS_index;  %A_set: active RRH set, D_set: inactive RRH set
        
        %% Process Deflation Procedure
        for A=0:L-1
            
            params.Inactive_index=D_set; 
            params.Active_index=A_set;
            params.weight=(1/eita)*ones(L,1);
            params.rankone=true; %return rankone solution
            [Wsolution, feasible] = powermin_cvx(params); %power minimization given the active RRH set
            
            if feasible==1
                
                TotalPower_IRLSAlternating_3(lq)=(1/eita)*norm(Wsolution,'fro')^2+sum(Pc(A_set));  %recode current values
                TransmitPower_IRLSAlternating_3(lq)=(1/eita)*norm(Wsolution,'fro')^2;
                A_number_IRLSAlternating_3(lq)=length(A_set);
                
                D_set=[D_set, BS_index(A+1)];  %%%process next selection: updata active RRH set and number
                A_set=setdiff([1:L],D_set); 
                
            else
                break;
            end
        end   
    end  
end
          
%% GSBF with Alternating Optimization Based Group Sparsity Design%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if GSBF_Alternating==true
    
    %% Compute Approximated Group Sparse Beamformer
    
    alternating_max=0; alternating_abs=10^99; 
    value1=10^99; value2=0;
    mu=1/L.*ones(L,1);
  
    while alternating_max<=20
        
        alternating_max=alternating_max+1; %interation numbers
        value1=value2; %recode the old objective value
        
        params.Inactive_index=[]; 
        params.Active_index=[1:L];
        params.weight=Pc./mu;
        params.rankone=false; %return general rank solution
        [Q_GSBF, feasible_GSBF] = powermin_cvx(params);
        
        if feasible_GSBF==false
            break;
        else
            value_temp=zeros(L,1);
            for l=1:L 
                temp1=0;
                for m=1:M
                    temp1=temp1+trace(Q_GSBF(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
                end
                value_temp(l)=real(temp1+10^(-3)/(M*sum(N_set)));
            end
            value_temp=(Pc./eita).*value_temp;
            value2=(1./mu)'*value_temp; %new objective value
            mu=sqrt(value_temp)./sum(sqrt(value_temp)); %updata mu                       
            alternating_abs=abs(value1-value2); %absolute value of the adjacent objective values
        end
    end
    
    if feasible_GSBF==false
        SU_channel_temp=0;  %This channel is infeasible
        break;
    else  
        %% RRH Ordering Based on the Approximated Group Sparse Beamformer
        Value_GSBF=zeros(L,1); 
        for l=1:L 
            temp1=0;
            for m=1:M
                temp1=temp1+trace(Q_GSBF(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
            end
            Value_GSBF(l)=temp1;
        end
        
        %% Sparsity Parttern Based RRH Ordering 
      %  [BS,BS_index]=sort(Value_GSBF);
        
        %% Proposed RRH Ordering
        ChannelGain=zeros(L,1);
        for l=1:L
            ChannelGain(l)=norm(vec(params.H(sum(N_set(1:(l-1)))+1:sum(N_set(1:l)),:,:)),'fro')^2;
        end
        [BS,BS_index]=sort(sqrt(ChannelGain.*eita./Pc).*Value_GSBF);
        
        %%
        D_set=[]; A_set=BS_index;  %A_set: active RRH set, D_set: inactive RRH set
        
        %% Process Deflation Procedure
        for A=0:L-1
            
            params.Inactive_index=D_set; 
            params.Active_index=A_set;
            params.weight=(1/eita)*ones(L,1);
            params.rankone=true; %return rankone solution
            [Wsolution, feasible] = powermin_cvx(params); %power minimization given the active RRH set
            
            if feasible==1
                
                TotalPower_GSBFAlternating(lq)=(1/eita)*norm(Wsolution,'fro')^2+sum(Pc(A_set));  %recode current values
                TransmitPower_GSBFAlternating(lq)=(1/eita)*norm(Wsolution,'fro')^2;
                A_number_GSBFAlternating(lq)=length(A_set);
                
                D_set=[D_set, BS_index(A+1)];  %%%process next selection: updata active RRH set and number
                A_set=setdiff([1:L],D_set); 
                
            else
                break;
            end
        end   
    end 
    
end


%% Baseline Algorith with l1/l_infty Group Sparsity Design%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Baseline==true
    
    %% Compute Approximated Group Sparse Beamformer
    params.Inactive_index=[]; 
    params.Active_index=[1:L];
    params.weight=(Pc(1)/eita*length(params.Active_index))*ones(length(params.Active_index),1);
    params.rankone=false; %return general rank solution
    [Q_GSBF, feasible_GSBF] = baseline_cvx(params);
    
    if feasible_GSBF==false
        SU_channel_temp=0;  %This channel is infeasible
        break;
    else  
        %% RRH Ordering Based on the Approximated Group Sparse Beamformer
        Value_GSBF=zeros(L,1); 
        for l=1:L 
            temp1=0;
            for m=1:M
                temp1=temp1+trace(Q_GSBF(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
            end
            Value_GSBF(l)=temp1;
        end
        
        %% Sparsity Parttern Based RRH Ordering 
        [BS,BS_index]=sort(Value_GSBF);

        
        D_set=[]; A_set=BS_index;  %A_set: active RRH set, D_set: inactive RRH set
        %% Process Deflation Procedure
        for A=0:L-1
            
            params.Inactive_index=D_set; 
            params.Active_index=A_set;
            params.weight=(1/eita)*ones(L,1);
            params.rankone=true; %return rankone solution
            [Wsolution, feasible] = powermin_cvx(params); %power minimization given the active RRH set
            
            if feasible==1
                
                TotalPower_Baseline(lq)=(1/eita)*norm(Wsolution,'fro')^2+sum(Pc(A_set));  %recode current values
                TransmitPower_Baseline(lq)=(1/eita)*norm(Wsolution,'fro')^2;
                A_number_Baseline(lq)=length(A_set);
                
                D_set=[D_set, BS_index(A+1)];  %%%process next selection: updata active RRH set and number
                A_set=setdiff([1:L],D_set); 
            else
                break;
            end
        end   
    end 
    
end

%% Exhaustive Search%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Exhaustive==true
    
    P_optimal=10^99; 
    D_set=[]; A_set=[1:L];  %A_set: active RRH set, D_set: inactive RRH set
    D_optimal=[]; A_optimal=[1:L];  % optimal active and inactive RRH set
    
    for A=0:L-1
        BS_select=nchoosek([1:L],A);
        feasible_temp=0; %recorder feasible times
        for s=1:nchoosek(L,A)
            
            A_set=setdiff([1:L],BS_select(s,:));   %A_set: active RRH set, D_set: inactive RRH set
            D_set=setdiff([1:L], A_set);
            
            params.Inactive_index=D_set; 
            params.Active_index=A_set;
            params.weight=(1/eita)*ones(L,1);
            params.rankone=true; %return rankone solution
            [Wsolution, feasible] = powermin_cvx(params); %power minimization given the active RRH set
            
            if  feasible==1   %%%if feasible
                temp1=(1/eita)*norm(Wsolution,'fro')^2+sum(Pc(A_set));  %network power consumption
                feasible_temp=feasible_temp+1;
                if P_optimal>temp1
                    P_optimal=temp1;   %update optimal network power, active RRH set
                    D_optimal=D_set; A_optimal=A_set;  %record optimal set
                end
            end
        end
        
        if feasible_temp==0  %all the available sets are infeasible
            break;
        end
    end
    
%% Output for the Particular QoS
if P_optimal<10^50
    TotalPower_Exhaustive(lq)=P_optimal;  %Network Power consumption
    TransmitPower_Exhaustive(lq)=P_optimal-sum(Pc(A_optimal)); %Transmit Power Consumption
    A_number_Exhaustive(lq)=length(A_optimal);  %Number of active RRHs
else
    SU_channel_temp=0;  % this channel if infeasible
    break;
end
end
end
    
%% Output for the Particular Channel Realization
if SU_channel_temp>0
    SU_counter=SU_counter+1; 
    if CBF==true
        TotalPower_CBF_temp=TotalPower_CBF_temp+TotalPower_CBF;
        TransmitPower_CBF_temp=TransmitPower_CBF_temp+TransmitPower_CBF;
    end
    if IRLS==true
        TotalPower_IRLSAlternating_temp=TotalPower_IRLSAlternating_temp+TotalPower_IRLSAlternating;
        TransmitPower_IRLSAlternating_temp=TransmitPower_IRLSAlternating_temp+TransmitPower_IRLSAlternating; 
        A_number_IRLSAlternating_temp=A_number_IRLSAlternating_temp+A_number_IRLSAlternating;
    end 
    if IRLS_2==true
        TotalPower_IRLSAlternating_2_temp=TotalPower_IRLSAlternating_2_temp+TotalPower_IRLSAlternating_2;
        TransmitPower_IRLSAlternating_2_temp=TransmitPower_IRLSAlternating_2_temp+TransmitPower_IRLSAlternating_2; 
        A_number_IRLSAlternating_2_temp=A_number_IRLSAlternating_2_temp+A_number_IRLSAlternating_2;
    end 
    if IRLS_3==true
        TotalPower_IRLSAlternating_3_temp=TotalPower_IRLSAlternating_3_temp+TotalPower_IRLSAlternating_3;
        TransmitPower_IRLSAlternating_3_temp=TransmitPower_IRLSAlternating_3_temp+TransmitPower_IRLSAlternating_3; 
        A_number_IRLSAlternating_3_temp=A_number_IRLSAlternating_3_temp+A_number_IRLSAlternating_3;
    end 
    if GSBF_Alternating==true
        TotalPower_GSBFAlternating_temp=TotalPower_GSBFAlternating_temp+TotalPower_GSBFAlternating;
        TransmitPower_GSBFAlternating_temp=TransmitPower_GSBFAlternating_temp+TransmitPower_GSBFAlternating; 
        A_number_GSBFAlternating_temp=A_number_GSBFAlternating_temp+A_number_GSBFAlternating;
    end  
    if Baseline==true
        TotalPower_Baseline_temp=TotalPower_Baseline_temp+TotalPower_Baseline;
        TransmitPower_Baseline_temp=TransmitPower_Baseline_temp+TransmitPower_Baseline; 
        A_number_Baseline_temp=A_number_Baseline_temp+A_number_Baseline;
    end
    
    if Exhaustive==true
        TotalPower_Exhaustive_temp=TotalPower_Exhaustive_temp+TotalPower_Exhaustive;
        TransmitPower_Exhaustive_temp=TransmitPower_Exhaustive_temp+TransmitPower_Exhaustive; 
        A_number_Exhaustive_temp=A_number_Exhaustive_temp+A_number_Exhaustive;
    end
    
end
    disp('channel complete');
end

%% Final Results
if CBF==true
    TotalPower_CBF=TotalPower_CBF_temp./SU_counter;
    TransmitPower_CBF=TransmitPower_CBF_temp./SU_counter;
end
if IRLS==true
    TotalPower_IRLSAlternating=TotalPower_IRLSAlternating_temp./SU_counter;
    TransmitPower_IRLSAlternating=TransmitPower_IRLSAlternating_temp./SU_counter;
    A_number_IRLSAlternating=A_number_IRLSAlternating_temp/SU_counter;
end
if IRLS_2==true
    TotalPower_IRLSAlternating_2=TotalPower_IRLSAlternating_2_temp./SU_counter;
    TransmitPower_IRLSAlternating_2=TransmitPower_IRLSAlternating_2_temp./SU_counter;
    A_number_IRLSAlternating_2=A_number_IRLSAlternating_2_temp/SU_counter;
end
if IRLS_3==true
    TotalPower_IRLSAlternating_3=TotalPower_IRLSAlternating_3_temp./SU_counter;
    TransmitPower_IRLSAlternating_3=TransmitPower_IRLSAlternating_3_temp./SU_counter;
    A_number_IRLSAlternating_3=A_number_IRLSAlternating_3_temp/SU_counter;
end
if GSBF_Alternating==true
    TotalPower_GSBFAlternating=TotalPower_GSBFAlternating_temp./SU_counter;
    TransmitPower_GSBFAlternating=TransmitPower_GSBFAlternating_temp./SU_counter;
    A_number_GSBFAlternating=A_number_GSBFAlternating_temp/SU_counter;
end

if Baseline==true
    TotalPower_Baseline=TotalPower_Baseline_temp./SU_counter;
    TransmitPower_Baseline=TransmitPower_Baseline_temp./SU_counter;
    A_number_Baseline=A_number_Baseline_temp/SU_counter;
end

if Exhaustive==true
    TotalPower_Exhaustive=TotalPower_Exhaustive_temp./SU_counter;
    TransmitPower_Exhaustive=TransmitPower_Exhaustive_temp./SU_counter;
    A_number_Exhaustive=A_number_Exhaustive_temp/SU_counter;
end
%% Plot everage network power consumption%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
if CBF==true
    plot(Q,TotalPower_CBF,'b-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end
if IRLS==true
    plot(Q,TotalPower_IRLSAlternating,'r-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end
if IRLS_2==true
    plot(Q,TotalPower_IRLSAlternating_2,'m-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end
if IRLS_3==true
    plot(Q,TotalPower_IRLSAlternating_3,'c-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end
if GSBF_Alternating==true
    plot(Q,TotalPower_GSBFAlternating,'y-d','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end

if Baseline==true
    plot(Q,TotalPower_Baseline,'black-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end

if Exhaustive==true
    plot(Q,TotalPower_Exhaustive,'g-d','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end

h=legend('CBF','IRLS for network minimization with p=1','IRLS for network minimization with p=0.5', 'Baseline','Exhaustive Search');  
xlabel('Target SINR [dB]','fontsize',12,'fontweight','b','fontname','helvetica');
ylabel('Average Power Consumption [W]','fontsize',14,'fontweight','b','fontname','helvetica');


%% Plot everage transmit power consumption%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
if CBF==true
    plot(Q,TransmitPower_CBF,'b-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end
if IRLS==true
    plot(Q,TransmitPower_IRLSAlternating,'r-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end
if IRLS_2==true
    plot(Q,TransmitPower_IRLSAlternating_2,'m-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end
if IRLS_3==true
    plot(Q,TransmitPower_IRLSAlternating_3,'c-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end
if GSBF_Alternating==true
    plot(Q,TransmitPower_GSBFAlternating,'y-d','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end

if Baseline==true
    plot(Q,TransmitPower_Baseline,'black-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end

if Exhaustive==true
    plot(Q,TransmitPower_Exhaustive,'g-d','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end

h=legend('CBF','IRLS for network minimization with p=1','IRLS for network minimization with p=0.5', 'Baseline','Exhaustive Search');  
xlabel('Target SINR [dB]','fontsize',12,'fontweight','b','fontname','helvetica');
ylabel('Average Transmit Power Consumption [W]','fontsize',14,'fontweight','b','fontname','helvetica');


%% Plot everage number of active RRHs%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
if CBF==true
    plot(Q,L*ones(length(Q),1),'b-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end
if IRLS==true
    plot(Q,A_number_IRLSAlternating,'r-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end
if IRLS_2==true
    plot(Q,A_number_IRLSAlternating_2,'m-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end
if IRLS_3==true
    plot(Q,A_number_IRLSAlternating_3,'c-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end
if GSBF_Alternating==true
    plot(Q,A_number_GSBFAlternating,'y-d','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end

if Baseline==true
    plot(Q,A_number_Baseline,'black-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end

if Exhaustive==true
    plot(Q,A_number_Exhaustive,'g-d','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end

h=legend('CBF','IRLS for network minimization with p=1','IRLS for network minimization with p=0.5', 'Baseline','Exhaustive Search');  

xlabel('Target SINR [dB]','fontsize',12,'fontweight','b','fontname','helvetica');
ylabel('Average Number of Active RRHs','fontsize',14,'fontweight','b','fontname','helvetica');
save('matlab1.mat');