clc;clear all;
%cvx_solver sedumi
%cvx_solver sdpt3
%cvx_solver mosek
cvx_solver sdpt3
%cvx_solver_settings('max_iters', 10^5, 'scale', 1,'eps',1.00e-5);
%cvx_solver_settings('SCALE', 1);
cvx_quiet(false)

%% Problem Data (I)
channel_generating=1;
LC=1; %# loops for channel realizarions
Area=10*10^3;
L=9; %Number of RAUs
N1=2;  % Antennas of each RAU

M=10; %Number of Multicast Groups
K1=2; %Number of Mobile Users in Each Multicast Group

%Q=[8:1:11];  %QoS in dB
Q=8;
N_set=N1*ones(L,1); % Set of Antennas for all the RAU
N=sum(N_set);
K_set=K1*ones(M,1); %Set of Numbers of Mobile Users in the Multicast Groups
K=sum(K_set);


Pc=3*ones(L,1)+[0:1:(L-1)]';  %power consumption of the fronthaul link
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
params.r_set=10^(Q/10)*ones(M,1);

if channel_generating==true
for ss=1:LC   %generate channel
    for m=1:length(K_set)
        D(:,m,:)=largefading(Area,L,K_set(m),N_set);
    end
    for m=1:1:length(K_set)
        H(:,m,:,ss)=smallfading(Area,L,K_set(m),N_set, D(:,m,:));
    end
end
save('H.mat','H');
end
load('H.mat');
params.alpha=1;
%% fesible detection
 params.H=H;
     
ebs=10^(-3);
params.alpha=1;
K_activeset=K_set;
M_activeindex=[1:M];
K_index=[1:1:sum(K_set)]; %%user can be supported
w=ones(sum(K_set),1);

iterative_max=0; iterative_abs=10^99;
value1=0;value2=10^99;

 params.K_activeset=K_activeset;
 params.M_activeindex=M_activeindex;
 params.K_index=K_index;
 params.rankone=false;
 params.weight=1/eita*ones(L,1);
 
 
 objective=[];
 index=1:1:21;
while iterative_max<=20
    iterative_max=iterative_max+1; %interation numbers
    value1=value2; %recode the old objective value
    params.w=w;
    
    [x,W,feasible]=powermin_uac(params);
    
    value_temp=zeros(L,1);
    for l=1:L 
        temp1=0;
        for m=1:M
            value_temp(l)=value_temp(l)+eita*trace(W(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
        end
    end
    value2=(1-params.alpha)*sum(value_temp)+params.alpha*w'*(abs(x).^2);
    w=1./(abs(x).^2+ebs);
    objective(iterative_max)=value2;
    iterative_abs=abs(value2-value1);
    
end


    plot(index,objective,'b-o','LineWidth',1.5, 'MarkerSize',10)
