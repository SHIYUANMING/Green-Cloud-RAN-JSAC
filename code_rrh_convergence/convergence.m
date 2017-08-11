clear all;
clc
% 
cvx_solver sdpt3;
% % % 
%cvx_solver_settings('MAX_ITERS', 10^4, 'SCALE',1);
channel_generation=0; 
 
%% Problem Data
Area =2*(10^3);
L=6;  %%Number of RAUs
NL=2; %%Number of Antennas of each RAU
N=L*NL;  %%Number of total antennas 
N_set=NL*ones(L,1); % Set of Antennas for all the RAU

M=2;  %%Number of Multicast Group
MK=2;  %%Number of users of each Multicast Group
K=M*MK;  %%Number of total users
K_set=MK*ones(M,1); %Set of Numbers of Mobile Users in the Multicast Groups
   
Pc=5.6*ones(L,1)+0*[0:1:(L-1)]';  %power consumption of the fronthaul link
eita=1/4;
%e=0.001;
e=10^(-3);
p=1;
%QoS=[0:2:10]';  %%QoS in dB
Q=4;
LC=1; %$ loops for channel realizations


params.L=L;
params.N=N;
params.M=M;
params.K=K;
params.K_set=K_set;
params.N_set=N_set;
params.delta_set=ones(M,1);
params.P_set=10^(0)*ones(L,1);
params.r_set=10^(Q/10)*ones(M,1);
%%%%%%%%%%%%%%%%%%%Generate Channel%%%%%%%%%%%%%%%%%
% for lc=1:1:LC
%     for m=1:1:length(K_set)
%         
%         H(:,m,:,lc)=channel_realization(Area,L,K_set(m),N_set);
%     end
% end

%% Wyner channel model
if channel_generation==1
  for lc=1:1:LC
    for m=1:length(K_set)
        D(:,m,:)=largefading(Area,L,K_set(m),N_set);
    end
    
    for m=1:1:length(K_set)
        H(:,m,:,lc)=smallfading(Area,L,K_set(m),N_set, D(:,m,:));
    end
 end
%save('H.mat','H');
end
load('H.mat');
params.H=H;  
   
alternating_max=31; alternating_abs=10^99;
alternating=0;
value1=10^99; value2=0;
w=ones(L,1);
objective=[];
index=1:1:alternating_max+1;
    while alternating<=alternating_max
        
        alternating=alternating+1; %interation numbers
        value1=value2; %recode the old objective value
        
        params.Inactive_index=[]; 
        params.Active_index=[1:L];
        params.w=w;
        params.rankone=false; %return general rank solution
        [Q_GSBF, feasible_GSBF] = powermin_cvx(params);
        
        if feasible_GSBF==false
            break;
        else
            value_temp=zeros(L,1);
            for l=1:L 
                for m=1:M
                    value_temp(l)=value_temp(l)+trace(Q_GSBF(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
                end                
            end
            
            objective(alternating)=sqrt(Pc(1)/eita)*sum((value_temp+e^2).^(p/2));
            
            w=Pc.*(value_temp+e^2).^(p/2-1);
            
            %%
            if e<=10^(-8)
                e=e/10;
            end
            %
            alternating_abs(alternating)=abs(value1-value2); %absolute value of the adjacent objective values
        end
    end
    %semilogy(index,alternating_abs,'b-o','LineWidth',1.5, 'MarkerSize',10)
    %figure;
     plot([0:alternating_max],objective./objective(1),'g-o','LineWidth',1.5, 'MarkerSize',10);
     hold on;
conver=objective./objective(1);
conver_data=zeros(length(conver),2);
conver_data(:,1)=[0:alternating_max];
conver_data(:,2)=conver;
save('conver_data.dat','conver_data','-ascii')

       