function [Vsolution,feasible] = powermin_cvx(params)

%prob_to_socp: maps PARAMS into a struct of SOCP matrices
%input struct 'parms' has the following fields:
%params.L;   %'L': # RRHs
%params.K;    %'K': # MUs
%params.N_set;  %set of antennas at all the RRHs
%params.Inactive_index;  %Index of the inactive RRHs
%params.Active_index;   %index of the active RRHs
%params.delta_set; %set of noise covariance

%%%%%%%%%%%%%%Problem Instances%%%%%%%%%%%%%
%params.r_set;  %set of SINR thresholds
%params.H;  %Channel Realization
%params.P_set;   %set of transmit power constraints at all the RRHs

%%%%%%%%Problem Data%%%%%%%
K_set=params.K_set;   %Mx1 vector: Set of Numbers of Mobile Users in M Multicast Group
r_set=params.r_set;   %Mx1 vector: QoS Requirements of Mobile Users for each Muliticast Group: assuming all the users have the same QoS requirments in the same group
delta_set=params.delta_set;  %Mx1 vector: noise covariance: assuming same value in the same group

Inactive_index=params.Inactive_index;  %Index of the inactive RRHs
Active_index=params.Active_index;       %Index of the inactive RRHs
N_set=params.N_set;
A_set=params.N_set(Active_index);   %Lx1 vector: active RAU antennas set
P_set=params.P_set(Active_index);  %Lx1 vector: active RAU transmit power set
Nl_set=[];%% total antennas of RAUs
Nl=0;
for l=1:1:length(Active_index)
    for ll=1:1:A_set(l)
        if Active_index(l)==1
            Nl_set=[Nl_set,ll];
        else
         Nl_set=[Nl_set,sum(N_set(1:Active_index(l)-1))+ll];
        end
    end
end
H=params.H(Nl_set,:,:);  %NxMxK channel matrix: N=sum(N_set), M=length(K_set), K equals the number of mobile users in each group (assuming each group has the same number of mobile users)

weight=params.weight(Active_index);  %length(Active_index)x1: weigth for each group beamformer
rankone=params.rankone;  %reture rankone solution or not

%%%%%%%%CVX+SCS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cvx_begin

variable Q(sum(A_set), sum(A_set), length(K_set)) hermitian semidefinite;   %Variable for N x N x M beamforming matrix for M groups
variable t(length(Active_index),1);  %slack variables for each beamforming matrix for the RAUs
minimize (weight'*t)  % group sparsity inducing minimization
subject to

%% RAUs Transmit Power Constraints
%% Active RAUs
for l=1:length(Active_index)  
    temp1=0;
    for m=1:length(K_set)
        temp1=temp1+trace(Q(sum(A_set(1:l-1))+1:sum(A_set(1:l)),sum(A_set(1:l-1))+1:sum(A_set(1:l)),m));
    end
    temp1<=P_set(l);
    temp1<=t(l);
    t(l)>=0;
end

%% Inactive RAUs
% for l=1:length(Inactive_index)  
%     temp1=0;
%     for m=1:length(K_set)
%         temp1=temp1+trace(Q(sum(N_set(1:Inactive_index(l)-1))+1:sum(N_set(1:Inactive_index(l))),sum(N_set(1:Inactive_index(l)-1))+1:sum(N_set(1:Inactive_index(l))),m));
%     end
%     temp1==0;
% end


%% QoS Constraints
Q_sum=0;
for m=1:length(K_set)
    Q_sum=Q_sum+Q(:,:,m);
end

for m=1:length(K_set)
    for k=1:K_set(m)
        G=(1+r_set(m))*Q(:,:,m)-r_set(m)*Q_sum;
        (H(:,m,k)'*G*H(:,m,k)-r_set(m)*delta_set(m)^2)==hermitian_semidefinite(1);
    end
end
cvx_end

 
%% build output
%%
     if  strfind(cvx_status,'Solved') 
         feasible=true;
        
         
         if rankone==true
         Vsolution=Rankone(Q);
         else
         Vsolution=Q;
         end
         
     else
         feasible=false;
         Vsolution=[];
     end