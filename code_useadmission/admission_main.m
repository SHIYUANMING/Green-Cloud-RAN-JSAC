clc;clear all;
%cvx_solver sedumi
%cvx_solver sdpt3
%cvx_solver mosek
cvx_solver sdpt3
%cvx_solver_settings('max_iters', 10^5, 'scale', 1,'eps',1.00e-5);
%cvx_solver_settings('SCALE', 1);
cvx_quiet(true)

%% Problem Data (I)
channel_generating=1;
LC=2; %# loops for channel realizarions


WEIGHTED_L2_P1=1;
WEIGHTED_L2_P2=1;
WEIGHTED_L2_P3=0;
WEIGHTED_L1=1;
L1=1;
MDR=1;
Exhaustive=1;

%%
Area=10*10^3;
L=6; %Number of RAUs
N1=2;  % Antennas of each RAU

M=4; %Number of Multicast Groups
K1=2; %Number of Mobile Users in Each Multicast Group

Q=[4:1:10];  %QoS in dB
%Q=8;

%%
N_set=N1*ones(L,1); % Set of Antennas for all the RAU
N=sum(N_set);
K_set=K1*ones(M,1); %Set of Numbers of Mobile Users in the Multicast Groups
%K_set=[1,1,1,1,2,2];
K=sum(K_set);
for i=1:1:M
    K_sum(i)=sum(K_set(1:i));
end

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

%%%%%%initialize results%%%%%%%%%%%%%%%%%%%%%%

user_temp_weightedl2_p1=0;%%our algorithm
user_temp_weightedl2_p2=0;
user_temp_weightedl2_p3=0;
user_temp_weightedl1=0;
user_temp_l1=0;
user_temp_mdr=0;%%
user_temp_exhaustive=0;
transmitpower_temp_weightedl2_p1=0;
transmitpower_temp_weightedl2_p2=0;
transmitpower_temp_weightedl2_p3=0;
transmitpower_temp_weightedl1=0;
transmitpower_temp_l1=0;
transmitpower_temp_mdr=0;
transmitpower_temp_exhaustive=0;
if channel_generating==true
    D=zeros(L,M,max(K_set));
    H=zeros(sum(N_set),M,max(K_set),LC);
for ss=1:LC   %generate channel
    for m=1:length(K_set)
        D(:,m,1:K_set(m))=largefading(Area,L,K_set(m),N_set);
    end
    for m=1:1:length(K_set)
        H(:,m,1:K_set(m),ss)=smallfading(Area,L,K_set(m),N_set, D(:,m,1:K_set(m)));
    end
end
save('H.mat','H');
%% standard channel
% for lc=1:1:LC
%     for m=1:1:length(K_set)
%         H(:,m,1:K_set(m),lc)=channel_realization(Area,L,K_set(m),N_set);
%     end
% end
% save('H.mat','H');
end
load('H.mat');
%% params for MDR
 params.epsilon=0;
 params.d=0.0001;
%% params for our algorithm
params.alpha=1;
p1=1;p2=0.5;p3=0;

SU_counter=0;

for lp=1:LC
      
    user_set_weightedl2_p1=zeros(length(Q),1);
    user_set_weightedl2_p2=zeros(length(Q),1);
    user_set_weightedl2_p3=zeros(length(Q),1);
    user_set_weightedl1=zeros(length(Q),1);
    user_set_l1=zeros(length(Q),1);
    user_set_mdr=zeros(length(Q),1);  
    user_set_exhaustive=zeros(length(Q),1);
    transmitpower_weightedl2_p1=zeros(length(Q),1);
    transmitpower_weightedl2_p2=zeros(length(Q),1);
    transmitpower_weightedl2_p3=zeros(length(Q),1);
    transmitpower_weightedl1=zeros(length(Q),1);
    transmitpower_l1=zeros(length(Q),1);
    transmitpower_mdr=zeros(length(Q),1);
    transmitpower_exhaustive=zeros(length(Q),1);
    params.H=H(:,:,:,lp);  %NxMxK channel matrix
    SU_flag=0;
for lq=1:length(Q)
    params.r_set=10^(Q(lq)/10)*ones(M,1);  %set of SINR thresholds  for each group 
    params.Active_index=[1:1:L];
    params.Inactive_index=[];
    params.weight=1/eita*ones(L,1);
    params.rankone=false;
    [Vsolution,Vfeasible]=powermin_cvx(params);  %%detect the whether the original algorithm is infeasible   
    if Vfeasible==1
        SU_flag=1;
        break
    end

%% User admission weighted l2norm with p=1
if WEIGHTED_L2_P1==1
ebs=10^(-3);

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
while iterative_max<=20
    iterative_max=iterative_max+1; %interation numbers
    value1=value2; %recode the old objective value
    params.w=w;
    
    [x,W,feasible]=powermin_uacl2(params);
    
    value_temp=zeros(L,1);
    for l=1:L 
        temp1=0;
        for m=1:M
            value_temp(l)=value_temp(l)+eita*trace(W(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
        end
    end
    value2=(1-params.alpha)*sum(value_temp)+params.alpha*w'*(abs(x).^2);
    if iterative_max<=6
    expindex=iterative_max;
    else expindex=5;
    end
    w=(abs(x).^2+ebs).^(p1/2-1);
    iterative_abs=abs(value2-value1);
end
%%user selection  
x_max=0;%%recoded during under new qos constriants
    for i=1:1:K  
        [x_max,sinry]=max(x);
        M_index_original=find(K_sum>=K_index(sinry));
        sinrx=find(M_activeindex==M_index_original(1));%%the index of multicast group after ignoring some.
        K_activeset(sinrx)=K_activeset(sinrx)-1;
        K_index=setdiff(K_index,K_index(sinry));
        x(sinry)=[];         
        if K_activeset(sinrx)==0
            K_activeset(sinrx)=[];
            M_activeindex(sinrx)=[];
        end%%%
        if sum(K_activeset)==0
            SU_flag=1;break;
        end
        params.K_activeset=K_activeset;
        params.M_activeindex=M_activeindex;
        params.K_index=K_index;
        params.rankone=false;
        params.weight=1/eita*ones(L,1);
        [Vsolution,feasible]=powermin_cvx_user(params);
        if feasible==1 
            user_set_weightedl2_p1(lq)=sum(K_activeset);
            V=Rankone(Vsolution);
            temp=0;
            for v=1:1:length(K_activeset)
                temp=temp+(1/eita)*norm(V(:,v))^2;
            end
            transmitpower_weightedl2_p1(lq)=temp;
            break;
        end         
    end
    if SU_flag==1
    break;
    end
end

%% User admission weighted l2 norm p=0.5
if WEIGHTED_L2_P2==1
ebs=10^(-3);

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
while iterative_max<=20
    iterative_max=iterative_max+1; %interation numbers
    value1=value2; %recode the old objective value
    params.w=w;
    
    [x,W,feasible]=powermin_uacl2(params);
    
    value_temp=zeros(L,1);
    for l=1:L 
        temp1=0;
        for m=1:M
            value_temp(l)=value_temp(l)+eita*trace(W(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
        end
    end
    value2=(1-params.alpha)*sum(value_temp)+params.alpha*w'*(abs(x).^2);
    if iterative_max<=6
    expindex=iterative_max;
    else expindex=5;
    end
    w=(abs(x).^2+ebs).^(p2/2-1);
    iterative_abs=abs(value2-value1);
end
%%user selection  
x_max=0;%%recoded during under new qos constriants
    for i=1:1:K  
        [x_max,sinry]=max(x);
        M_index_original=find(K_sum>=K_index(sinry));
        sinrx=find(M_activeindex==M_index_original(1));%%the index of multicast group after ignoring some.
        K_activeset(sinrx)=K_activeset(sinrx)-1;
        K_index=setdiff(K_index,K_index(sinry));
        x(sinry)=[];         
        if K_activeset(sinrx)==0
            K_activeset(sinrx)=[];
            M_activeindex(sinrx)=[];
        end%%%
        if sum(K_activeset)==0
            SU_flag=1;break;
        end
        params.K_activeset=K_activeset;
        params.M_activeindex=M_activeindex;
        params.K_index=K_index;
        params.rankone=false;
        params.weight=1/eita*ones(L,1);
        [Vsolution,feasible]=powermin_cvx_user(params);
        if feasible==1 
            user_set_weightedl2_p2(lq)=sum(K_activeset);
             V=Rankone(Vsolution);
            temp=0;
            for v=1:1:length(K_activeset)
                temp=temp+(1/eita)*norm(V(:,v))^2;
            end
            transmitpower_weightedl2_p2(lq)=temp;
            break;
        end         
    end
    if SU_flag==1
    break;
    end
end


%% User admission weighted l2 norm p=0
if WEIGHTED_L2_P3==1
ebs=10^(-3);

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
while iterative_max<=20
    iterative_max=iterative_max+1; %interation numbers
    value1=value2; %recode the old objective value
    params.w=w;
    
    [x,W,feasible]=powermin_uacl2(params);
    
    value_temp=zeros(L,1);
    for l=1:L 
        temp1=0;
        for m=1:M
            value_temp(l)=value_temp(l)+eita*trace(W(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
        end
    end
    value2=(1-params.alpha)*sum(value_temp)+params.alpha*w'*(abs(x).^2);
    if iterative_max<=6
    expindex=iterative_max;
    else expindex=5;
    end
    w=(abs(x).^2+ebs).^(p3/2-1);
    iterative_abs=abs(value2-value1);
end
%%user selection  
x_max=0;%%recoded during under new qos constriants
    for i=1:1:K  
        [x_max,sinry]=max(x);
        M_index_original=find(K_sum>=K_index(sinry));
        sinrx=find(M_activeindex==M_index_original(1));%%the index of multicast group after ignoring some.
        K_activeset(sinrx)=K_activeset(sinrx)-1;
        K_index=setdiff(K_index,K_index(sinry));
        x(sinry)=[];         
        if K_activeset(sinrx)==0
            K_activeset(sinrx)=[];
            M_activeindex(sinrx)=[];
        end%%%
        if sum(K_activeset)==0
            SU_flag=1;break;
        end
        params.K_activeset=K_activeset;
        params.M_activeindex=M_activeindex;
        params.K_index=K_index;
        params.rankone=false;
        params.weight=1/eita*ones(L,1);
        [Vsolution,feasible]=powermin_cvx_user(params);
        if feasible==1 
            user_set_weightedl2_p3(lq)=sum(K_activeset);
            V=Rankone(Vsolution);
            temp=0;
            for v=1:1:length(K_activeset)
                temp=temp+(1/eita)*norm(V(:,v))^2;
            end
            transmitpower_weightedl2_p3(lq)=temp;
            break;
        end         
    end
    if SU_flag==1
    break;
    end
end
%% User admission weighted l1 norm
if WEIGHTED_L1==1
ebs=10^(-3);

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
while iterative_max<=20
    iterative_max=iterative_max+1; %interation numbers
    value1=value2; %recode the old objective value
    params.w=w;
    
    [x,W,feasible]=powermin_uacl1(params);
    
    value_temp=zeros(L,1);
    for l=1:L 
        temp1=0;
        for m=1:M
            value_temp(l)=value_temp(l)+eita*trace(W(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
        end
    end
    value2=(1-params.alpha)*sum(value_temp)+params.alpha*w'*abs(x);
    w=(abs(x)+ebs).^(-1);
    iterative_abs=abs(value2-value1);
end
%%user selection  
x_max=0;%%recoded during under new qos constriants
    for i=1:1:K  
        [x_max,sinry]=max(x);
        M_index_original=find(K_sum>=K_index(sinry));
        sinrx=find(M_activeindex==M_index_original(1));%%the index of multicast group after ignoring some.
        K_activeset(sinrx)=K_activeset(sinrx)-1;
        K_index=setdiff(K_index,K_index(sinry));
        x(sinry)=[];         
        if K_activeset(sinrx)==0
            K_activeset(sinrx)=[];
            M_activeindex(sinrx)=[];
        end%%%
        if sum(K_activeset)==0
            SU_flag=1;break;
        end
        params.K_activeset=K_activeset;
        params.M_activeindex=M_activeindex;
        params.K_index=K_index;
        params.rankone=false;
        params.weight=1/eita*ones(L,1);
        [Vsolution,feasible]=powermin_cvx_user(params);
        if feasible==1
            user_set_weightedl1(lq)=sum(K_activeset);
            V=Rankone(Vsolution);
            temp=0;
            for v=1:1:length(K_activeset)
                temp=temp+(1/eita)*norm(V(:,v))^2;
            end
            transmitpower_weightedl1(lq)=temp;
            break;
        end         
    end
    if SU_flag==1
    break;
    end
end
%% User admission with l1-norm 
if L1==1
ebs=10^(-3);

K_activeset=K_set;
M_activeindex=[1:M];
K_index=[1:1:sum(K_set)]; %%user can be supported
w=ones(sum(K_set),1);

 params.K_activeset=K_activeset;
 params.M_activeindex=M_activeindex;
 params.K_index=K_index;
 params.rankone=false;
 params.weight=1/eita*ones(L,1);
 params.w=w;
[x,W,feasible]=powermin_uacl1(params);
%%user selection  
x_max=0;%%recoded during under new qos constriants
    for i=1:1:K  
        [x_max,sinry]=max(x);
        M_index_original=find(K_sum>=K_index(sinry));
        sinrx=find(M_activeindex==M_index_original(1));%%the index of multicast group after ignoring some.
        K_activeset(sinrx)=K_activeset(sinrx)-1;
        K_index=setdiff(K_index,K_index(sinry));
        x(sinry)=[];         
        if K_activeset(sinrx)==0
            K_activeset(sinrx)=[];
            M_activeindex(sinrx)=[];
        end%%%
        if sum(K_activeset)==0
            SU_flag=1;break;
        end
        params.K_activeset=K_activeset;
        params.M_activeindex=M_activeindex;
        params.K_index=K_index;
        params.rankone=false;
        params.weight=1/eita*ones(L,1);
        [Vsolution,feasible]=powermin_cvx_user(params);
        if feasible==1 
            user_set_l1(lq)=sum(K_activeset);
            V=Rankone(Vsolution);
            temp=0;
            for v=1:1:length(K_activeset)
                temp=temp+(1/eita)*norm(V(:,v))^2;
            end
            transmitpower_l1(lq)=temp;
            break;
        end         
    end
    if SU_flag==1
    break;
    end
end


%% User admission MDR %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if MDR==1
sinr_min=0;%%recoded during under new qos constriants
K_activeset=K_set;%%active set of users in different multicast group
M_activeindex=[1:M];%%index of the multicast group which is not totally ignored
K_index=[1:1:sum(K_set)]; %%user index can be supported
w=ones(sum(K_set),1);
while sum(K_activeset)~=0
    params.K_activeset=K_activeset;
    params.M_activeindex=M_activeindex;
    params.K_index=K_index;
    params.rankone=true;
    params.weight=1/eita*ones(L,1);
    [W,feasible]=powermin_mdr(params);
    sinr=ones(length(K_activeset),sum(K_activeset))*10^(Q(lq)/10);
    for i=1:1:length(K_activeset)
        m=M_activeindex(i);
        for k=1:1:K_activeset(i)
            kk=(sum(K_activeset(1:i-1))+k);
            ki=K_index(sum(K_activeset(1:i-1))+k)-sum(K_set(1:m-1));
            inter_i=(norm(W'*H(:,m,ki,lp)))^2-(norm(W(:,i)'*H(:,m,ki,lp)))^2;
            if (norm(W(:,i)'*H(:,m,ki,lp)))^2/(inter_i+params.delta_set(m)^2)>=params.r_set(m)
                sinr(i,kk)=params.r_set(m);
            else sinr(i,kk)=(norm(W(:,i)'*H(:,m,ki,lp)))^2/(inter_i+params.delta_set(m)^2);
            end
        end
    end
    sinr_min=min(sinr(:));
    [sinrx,sinry]=find(sinr==sinr_min);
    if sinr_min~=10^(Q(lq)/10)
        K_activeset(sinrx)=K_activeset(sinrx)-1;
        K_index=setdiff(params.K_index,K_index(sinry));
        if K_activeset(sinrx)==0
            K_activeset(sinrx)=[];
            M_activeindex(sinrx)=[];
        end%%%
        if sum(K_activeset)==0
            SU_flag=1;break;
        end
    else user_set_mdr(lq)=sum(K_activeset);  
         
         temp=0;
         for v=1:1:length(K_activeset)
             temp=temp+(1/eita)*norm(W(:,v))^2;
         end
         transmitpower_mdr(lq)=temp;
        break;
            
    end
end

if SU_flag==1
    break;
end

end


%% Exhaustive Search%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Exhaustive==true
   params.Inactive_index=[]; 
   params.Active_index=[1:L]; 
   params.rankone=false;
   params.weight=1/eita*ones(L,1);
   Feasible_temp=0;
    for A=1:K  %%number of users we didnt support
        USER_select=nchoosek([1:K],A);
        for s=1:nchoosek(K,A) %%different sets
            
             K_activeset=K_set;%%every time we reset the vectors
             M_activeindex=[1:M];
             K_index=[1:1:sum(K_set)];
             
            for i=1:1:A %%different users in each set
                sinry=USER_select(s,i); %% index of the user
                M_index_original=find(K_sum>=K_index(sinry));
                sinrx=find(M_activeindex==M_index_original(1)); %% index of the multicast the specific user belongs to 
                K_activeset(sinrx)=K_activeset(sinrx)-1;
                if K_activeset(sinrx)==0
                K_activeset(sinrx)=[];
                M_activeindex(sinrx)=[];
                end
            end
            K_index=setdiff([1:sum(K_set)],USER_select(s,:));   
            
            
            params.K_index=K_index;
            params.K_activeset=K_activeset;
            params.M_activeindex=M_activeindex;
            if sum(K_activeset)==0
            SU_flag=1;break;
            end
            
            [Vsolution,feasible]=powermin_cvx_user(params);
            
            if  feasible==1   %%%if feasible
                user_set_exhaustive(lq)=sum(K_activeset);
                Feasible_temp=1;
                V=Rankone(Vsolution);
                temp=0;
                for v=1:1:length(K_activeset)
                    temp=temp+(1/eita)*norm(V(:,v))^2;
                end
                transmitpower_exhaustive(lq)=temp;
            break;
            end
        end
        if SU_flag==1
            break;
        end
        if Feasible_temp==1
            break;
        end 
         
    end
    if SU_flag==1
            break;
    end 
        
end


end





if SU_flag==0
 SU_counter=SU_counter+1;
user_temp_weightedl2_p1=user_temp_weightedl2_p1+user_set_weightedl2_p1;
user_temp_weightedl2_p2=user_temp_weightedl2_p2+user_set_weightedl2_p2;
user_temp_weightedl2_p3=user_temp_weightedl2_p3+user_set_weightedl2_p3;
user_temp_weightedl1=user_temp_weightedl1+user_set_weightedl1;
user_temp_l1=user_temp_l1+user_set_l1;
user_temp_mdr=user_temp_mdr+user_set_mdr;
user_temp_exhaustive=user_temp_exhaustive+user_set_exhaustive;


transmitpower_temp_weightedl2_p1=transmitpower_temp_weightedl2_p1+transmitpower_weightedl2_p1;
transmitpower_temp_weightedl2_p2=transmitpower_temp_weightedl2_p2+transmitpower_weightedl2_p2;
transmitpower_temp_weightedl2_p3=transmitpower_temp_weightedl2_p3+transmitpower_weightedl2_p3;
transmitpower_temp_weightedl1=transmitpower_temp_weightedl1+transmitpower_weightedl1;
transmitpower_temp_l1=transmitpower_temp_l1+transmitpower_l1;
transmitpower_temp_mdr=transmitpower_temp_mdr+transmitpower_mdr;
transmitpower_temp_exhaustive=transmitpower_temp_exhaustive+transmitpower_exhaustive;
end


disp('channel complete');
end

%% Final Results

user_set_weightedl2_p1=user_temp_weightedl2_p1/SU_counter;
user_set_weightedl2_p2=user_temp_weightedl2_p2/SU_counter;
user_set_weightedl2_p3=user_temp_weightedl2_p3/SU_counter;
user_set_weightedl1=user_temp_weightedl1/SU_counter;
user_set_l1=user_temp_l1/SU_counter;
user_set_mdr=user_temp_mdr/SU_counter;
user_set_exhaustive=user_temp_exhaustive/SU_counter;


transmitpower_weightedl2_p1=transmitpower_temp_weightedl2_p1/SU_counter;
transmitpower_weightedl2_p2=transmitpower_temp_weightedl2_p2/SU_counter;
transmitpower_weightedl2_p3=transmitpower_temp_weightedl2_p3/SU_counter;
transmitpower_weightedl1=transmitpower_temp_weightedl1/SU_counter;
transmitpower_l1=transmitpower_temp_l1/SU_counter;
transmitpower_mdr=transmitpower_temp_mdr/SU_counter;
transmitpower_exhaustive=transmitpower_temp_exhaustive/SU_counter;
%% Plot everage network power consumption%%%%%%%%%%%%
figure;

if WEIGHTED_L2_P1==1
plot(Q,user_set_weightedl2_p1,'r-d','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
hold on
end

if WEIGHTED_L2_P2==1
plot(Q,user_set_weightedl2_p2,'b-x','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
hold on
end
if WEIGHTED_L2_P3==1
plot(Q,user_set_weightedl2_p3,'y-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
hold on
end
if WEIGHTED_L1==1
plot(Q,user_set_weightedl1,'m-*','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
hold on
end
 if L1==1
plot(Q,user_set_l1,'g-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
hold on
 end
if MDR==1
plot(Q,user_set_mdr,'c-v','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
hold on;
end
if Exhaustive==1
plot(Q,user_set_exhaustive,'k-+','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
hold on;    
end
xlabel('Target SINR [dB]','fontsize',12,'fontweight','b','fontname','helvetica');
ylabel('Average Number of Users','fontsize',14,'fontweight','b','fontname','helvetica');
%legend('Weighted l2 norm with p=1','Weighted l1 norm','l1 Norm','MDR')
legend('Weighted l2 norm with p=1','Weighted l2 norm with p=0.5','Weighted l1 norm','l1 Norm','MDR','Exhaustive Search')
%legend('Weighted L1','l1 Norm','MDR','Exhaustive Search')

%% TransmitPower
figure;

if WEIGHTED_L2_P1==1
plot(Q,transmitpower_weightedl2_p1,'r-d','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
hold on
end

if WEIGHTED_L2_P2==1
plot(Q,transmitpower_weightedl2_p2,'b-x','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
hold on
end
if WEIGHTED_L2_P3==1
plot(Q,transmitpower_weightedl2_p3,'y-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
hold on
end
if WEIGHTED_L1==1
plot(Q,transmitpower_weightedl1,'m-*','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
hold on
end
 if L1==1
plot(Q,transmitpower_l1,'g-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
hold on
 end
if MDR==1
plot(Q,transmitpower_mdr,'c-v','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
hold on;
end
if Exhaustive==1
plot(Q,transmitpower_exhaustive,'k-+','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
hold on;    
end
xlabel('Target SINR [dB]','fontsize',12,'fontweight','b','fontname','helvetica');
ylabel('Average Transmit Power Consumption','fontsize',14,'fontweight','b','fontname','helvetica');
%legend('Weighted l2 norm with p=1','Weighted l1 norm','l1 Norm','MDR')
legend('Weighted l2 norm with p=1','Weighted l2 norm with p=0.5','Weighted l1 norm','l1 Norm','MDR','Exhaustive Search')
save('matlab.mat')