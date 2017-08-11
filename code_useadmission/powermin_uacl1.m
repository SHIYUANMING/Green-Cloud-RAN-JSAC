function [xsolution,Qsolution,feasible] = powermin_uacl1(params)

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
M_activeindex=params.M_activeindex; %%multicast组全部关闭的index
K_set=params.K_set;   %Mx1 vector: Set of Numbers of Mobile Users in M Multicast Group
r_set=params.r_set(M_activeindex);   %Mx1 vector: QoS Requirements of Mobile Users for each Muliticast Group: assuming all the users have the same QoS requirments in the same group
delta_set=params.delta_set(M_activeindex);  %Mx1 vector: noise covariance: assuming same value in the same group
K_activeset=params.K_activeset;%%更新的K_set，某些可能是1；


N_set=params.N_set;  %Lx1 vector: active RAU antennas set
P_set=params.P_set;  %Lx1 vector: active RAU transmit power set
K_index=params.K_index; %%the index of users can be supported

    
H=params.H(:,:,:);  %NxMxK channel matrix: N=sum(N_set), M=length(K_set), K equals the number of mobile users in each group (assuming each group has the same number of mobile users)
alpha=params.alpha;
w=params.w;
weight=params.weight;  %length(Active_index)x1: weigth for each group beamformer
rankone=params.rankone;  %reture rankone solution or not

%%%%%%%%CVX+SCS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cvx_begin

variable Q(sum(N_set), sum(N_set), length(K_activeset)) hermitian semidefinite;   %Variable for N x N x M beamforming matrix for M groups
variable x(sum(K_activeset),1);  %%slack variables for objective function
variable t(length(N_set),1);  %slack variables for each beamforming matrix for the RAUs
variable s(sum(K_activeset),1);  %%slack variables for the slack objective variables.
minimize ((1-alpha)*(weight'*t)+alpha*w'*s)  % group sparsity inducing minimization
%minimize (w'*s)
subject to

%% RAUs Transmit Power Constraints
%% Active RAUs
for l=1:length(N_set)  
    temp1=0;
    for m=1:length(K_activeset)
        temp1=temp1+trace(Q(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
    end
    temp1<=P_set(l);
    temp1<=t(l);
    t(l)>=0;
end

%% Slack Variables
for i=1:1:sum(K_activeset)
    temp=x(i);
    abs(temp)<=s(i);
end



%% QoS Constraints
Q_sum=0;
for m=1:length(K_activeset)
    Q_sum=Q_sum+Q(:,:,m);
end

for i=1:1:length(K_activeset)  
    m=M_activeindex(i);
    for k=1:K_activeset(i) 
        G=(1+r_set(i))*Q(:,:,i)-r_set(i)*Q_sum;
        (H(:,m,K_index(sum(K_activeset(1:i-1))+k)-sum(K_set(1:m-1)))'*G*H(:,m,K_index(sum(K_activeset(1:i-1))+k)-sum(K_set(1:m-1)))-r_set(i)*delta_set(i)^2+x(sum(K_activeset(1:i-1))+k))==hermitian_semidefinite(1);
    end
end
cvx_end

 
%% build output
%%
     if  strfind(cvx_status,'Solved') 
         feasible=true;
        
         xsolution=x;
         if rankone==true
         Qsolution=Rankone(Q);
         else
         Qsolution=Q;
         end
         
     else
         feasible=false;
         Qsolution=[];
         xsolution=[];
     end
     
     %%((abs(x)).^2)