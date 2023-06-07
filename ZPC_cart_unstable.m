% ZPC: run ZPC with model and without the model
% i.e., run 
% 1-Robust data driven Predictive control scheme (ZPC) 
% 2- Same ZPC while knowing the model (RMPC-zono)
%
% Inputs:
%    none
%
% Outputs:
%    saved workspace
%
% Example: 
%
% See also: ---

% Author:       Amr Alanwar, Yvonne Stürz 
% Written:      25-March-2021 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
addpath(genpath('..\tbxmanager\'))
addpath(genpath('contSet\'))
addpath(genpath('matrixSet\'))
addpath(genpath('C:\Program Files\Mosek\10.0\'))

rand('seed',4500);

clear all
close all

% dimension of x
dim_x = 4;
dim_u = 1;

% Nominal system
mp = 2;
mc = 10;
l = 1;
g= 9.81;

% define discrete system
dt = 0.1;

% System in cont time
A = [0, 0, 1, 0; 0, 0, 0, 1; 0, mp*g/mc, 0, 0; 0, -(mp+mc)*g/(mc*l), 0, 0];
% A = eye(4) + [0, 0, dt, 0; 0, 0, 0, dt; 0, dt*mp*g/mc, 0, 0; 0, -dt*(mp+mc)*g/(mc*l), 0, 0];
B_ss = [0; 0; 1/mc; 1/(mc*l)];
C = eye(dim_x);
D = zeros(dim_x, 1);

sys_c = ss(A, B_ss, C, D);
sys_d = c2d(sys_c, dt);
%number of trajectories
initpoints = 50;
%number of steps for each trajectory
Tf = 30;
steps = Tf/dt;
%Total number of samples
totalsamples = initpoints*steps;
%% initial set and input
%reference input
uref = 0;
%reference output (steady state states) 
% ref = inv(eye(dim_x)-sys_d.A)*sys_d.B*uref; 
% Hovering condition (inverse doesn't work because A is not invertible) 
ref = [0; 0; 0; 0];

%output constraint
% y_lb = [-10;2;-10;-10;-10]; 
% y_ub = [10;10;10;10;10]; 
y_lb = [-1;-pi/4;-0.1;-.1]; 
y_ub = -1.*y_lb;
% y_ub = [1;0.1;.1;0.1]; 
intc = interval(y_lb,y_ub);


%initial point
% y0 = zeros(dim_x,1);
y0 = [0.5; 0.05; 0; 0];

%initial zonotope to generate data
% X0 = zonotope([y0,25*diag(ones(dim_x,1))]);
X0 = zonotope([y0,zeros(dim_x,1)]);

%input zonotope
u_lb = 5;
U = zonotope([uref,u_lb]);

%noise zontope W (modeling noise)
%less noise
wfac=0.001;
%more noise
%wfac=0.1;
W = zonotope([zeros(dim_x,1),wfac*ones(dim_x,1)]);

for i=1:size(W.generators,2)
    vec=W.Z(:,i+1);
    GW{i}= [ vec,zeros(dim_x,totalsamples-1)];
    for j=1:totalsamples-1
        GW{j+i}= [GW{i+j-1}(:,2:end) GW{i+j-1}(:,1)];
    end
end
% matrix zonotpe of noise w (M_w)
Wmatzono= matZonotope(zeros(dim_x,totalsamples),GW);

%measurement noise
%less measurement noise
vfac = 0.002;
%more measurement noise
%vfac = 0.02;
V = zonotope([zeros(dim_x,1),vfac*ones(dim_x,1)]);
CV = zeros(dim_x,totalsamples);
for i=1:size(V.generators,2)
    vec=V.Z(:,i+1);
    GV{i}= [ vec,zeros(dim_x,totalsamples-1)];
    for j=1:totalsamples-1
        GV{j+i}= [GV{i+j-1}(:,2:end) GV{i+j-1}(:,1)];
    end
end
% matrix zonotpe of noise v (M_v)
Vmatzono= matZonotope(CV,GV);

AV = sys_d.A*V;
% matrix zonotpe of  Av (M_Av)
VAmatzono = sys_d.A*Vmatzono;

% Control at each step
u = zeros(dim_u, totalsamples);
% randomly choose constant inputs for each step / sampling time
for i=1:totalsamples
    u(:, i) = randPoint(U);
end

%generate data from different trajectories with noise
x0 = X0.center;
x(:,1) = x0;
index=1;

utraj = zeros(initpoints*dim_x, steps, dim_u);
% Show trajectories
% TODO: Change to nonlinear dynamics 
for j=1:dim_x:initpoints*dim_x
    x(j:j+dim_x-1,1) = randPoint(X0);
    x_v(j:j+dim_x-1,1) =  x(j:j+dim_x-1,1) + randPoint(V);
    
    for i=1:steps
        utraj(j,i,:) = u(:,index);
        x(j:j+dim_x-1,i+1) = sys_d.A*x(j:j+dim_x-1,i) + sys_d.B*u(:,index) + randPoint(W);
        x_v(j:j+dim_x-1,i+1) =  x(j:j+dim_x-1,i+1) + randPoint(V);
        index=index+1;
    end
end



%prepeare Y_+ Y_-
index_0 =1;
index_1 =1;
% u_mean_vec_0 = zeros(dim_u, initpoints*dim_x*steps, dim_u);
for j=1:dim_x:initpoints*dim_x
    for i=2:steps+1
        x_meas_vec_1_v(:,index_1) = x_v(j:j+dim_x-1,i);
        x_meas_vec_1(:,index_1) = x(j:j+dim_x-1,i);
        index_1 = index_1 + 1;
    end
    for i=1:steps
        % u_mean_vec_0(:,index_0, :) = utraj(j,i,:);
        u_mean_vec_0(:,index_0) = utraj(j,i,:);
        x_meas_vec_0(:,index_0) = x(j:j+dim_x-1,i);
        x_meas_vec_0_v(:,index_0) = x_v(j:j+dim_x-1,i);
        index_0 = index_0 + 1;
    end
end
% U_data is U_-, Y_0T is Y_- , Y_1T is Y_+
U_data = u_mean_vec_0(:,1:totalsamples); %same as u
Y_0T = x_meas_vec_0_v(:,1:totalsamples);
Y_1T = x_meas_vec_1_v(:,1:totalsamples);



t_vec = dt*(0:steps);
% plot simulated trajectory
figure;
subplot(2,1,1)
hold on; box on; plot(t_vec, x(1,:),'b'); xlabel('t'); ylabel('x'); axis equal;
subplot(2,1,2)
hold on; box on; plot(t_vec, x(2,:),'b'); xlabel('t'); ylabel('\theta'); axis equal;
figure;
subplot(2,1,1)
hold on; box on; plot(t_vec, x(3,:),'b'); xlabel('t'); ylabel('\dot{x}','interpreter','latex'); axis equal;
subplot(2,1,2)
hold on; box on; plot(t_vec, x(4,:),'b'); xlabel('t'); ylabel('\dot{\theta}','interpreter','latex'); axis equal;
close all;
disp(max(x(2:4:end)))

%prepare M_Sigma which is a set of [A B]
AB = (Y_1T + -1* Vmatzono + -1*Wmatzono+VAmatzono)*pinv([Y_0T;U_data]);


%double check if the true A B is part of M_Sigma
intAB11 = intervalMatrix(AB);
intAB1 = intAB11.int;
intAB1.sup >= [sys_d.A,sys_d.B]
intAB1.inf <= [sys_d.A,sys_d.B]
% check the rank of the data
rank = rank([Y_0T;U_data])




%% Compute ZPC problem
%Horizon N for ZPC
N = 2;
%define output cost matrix
Qy = 1e3*eye(dim_x); 
%control cost matrix
Qu = 0.001*eye(dim_u);


execTimeZPC=[];
execTimeRMPC=[];
% ZPC number of time steps
maxsteps = Tf/dt;
% chosen time step for plotting 
chosedtimestep = 10;
for timesteps = 1:maxsteps
    if timesteps == 1
        % set the initial output to y0
        y_t(:,timesteps) = y0;
        y_t_model(:,timesteps) = y0;
        YPred(:,1) = y0;
    end
    
    
    % sdpvar variables
    u = sdpvar(dim_u*ones(1,N),ones(1,N));
    y = sdpvar(dim_x*ones(1,N+1),ones(1,N+1));
    alpha_u = sdpvar(dim_u,N);
    sinf = sdpvar(dim_x*ones(1,N+1),ones(1,N+1));
    ssup = sdpvar(dim_x*ones(1,N+1),ones(1,N+1));
    R={};
    R{1} = zonotope([y_t(:,timesteps)]);
    %set the first constraint as y_t = current y
    Constraints = [y_t(:,timesteps) == y{1}];%,...
   
    
    for i = 1:N
        %compute the reachable set for ZPC
        card_cen = [R{i}.center;u{i}];
        genleni = size(R{i}.generators,2);
        card_zono = zonotope([card_cen,[R{i}.generators;zeros(dim_u,genleni)]]);
        ABcard = intervalMatrix(AB)* card_zono;
        R{i+1} = zonotope([ABcard.center,[ABcard.generators,W.generators,V.generators,AV.generators]]);%AB * card_zono + W_sdp;
        %convert R to interval
        %extract center
        c = R{i+1}.Z(:,1);
        
        %determine left and right limit of the reahable set (convert to
        %interval)
        delta = sum(abs(R{i+1}.Z),2) - abs(c);
        leftLimit{i} = c - delta;
        rightLimit{i} = c + delta;
        
        %specify the constraints
        Constraints = [Constraints,...
            u{i} == U.center + U.generators * alpha_u(:,i) ,...
            y{i+1} - sinf{i} == leftLimit{i},...
            y{i+1} + ssup{i} == rightLimit{i},...
            y{i+1}   - sinf{i} >= intc.inf,...
            y{i+1}   + ssup{i} <= intc.sup,...
            sinf{i} >= zeros(dim_x,1),...
            ssup{i} >= zeros(dim_x,1),...
            alpha_u(:,i) <= ones(dim_u, 1) , ...
            alpha_u(:,i) >= -1*ones(dim_u, 1), ...
            ];
    end
    
    % chose the cost of ZPC
    Cost=0;
    for i=1:N
        Cost = Cost + (y{i+1}-ref)'*Qy*(y{i+1}-ref)+ (u{i}-uref)'*Qu*(u{i}-uref);
    end
    %solve ZPC
    options = sdpsettings('verbose',0,'solver','mosek');
    tic
    Problem = optimize(Constraints,Cost,options);
    execTimeZPC=[execTimeZPC,toc];
    Objective = double(Cost);
    uPred(:, timesteps) = double(u{1});
    YPred(:,timesteps+1) = double(y{2});
    %%

    %% save for plotting
    Rplotall{timesteps}= interval(zonotope([ double(R{2}.center), double(R{2}.generators)]));
    %%  ploting
    if chosedtimestep == timesteps
        for i =1:N+1
            RoverN{i}= zonotope([ double(R{i}.center), double(R{i}.generators)]) ;
            RoverN_int{i} = interval(RoverN{i});
            yoverN{i} =double(y{i});
            if i<N+1
                uoverN{i} =double(u{i});
            end
        end
    end
	
	
	
    %% ZPC given the model (RMPC-zono)
    % Control					
    alpha_u = sdpvar(dim_u,N);
    sinf = sdpvar(dim_x*ones(1,N+1),ones(1,N+1));
    ssup = sdpvar(dim_x*ones(1,N+1),ones(1,N+1));
    R={};
    R{1} = zonotope([y_t_model(:,timesteps)]);
    u_model = sdpvar(dim_u*ones(1,N),ones(1,N));
    y_model = sdpvar(dim_x*ones(1,N+1),ones(1,N+1));
    Constraints = [y_t_model(:,timesteps) == y_model{1}];
    for i = 1:N
        %card_cen = [y{i};u_model{i}];
        card_cen = [R{i}.center;u_model{i}];
        genleni = size(R{i}.generators,2);
        card_zono = zonotope([card_cen,[R{i}.generators;zeros(dim_u,genleni)]]);
        % give it true A B
        ABcard = [sys_d.A , sys_d.B]* card_zono;
        R{i+1} = zonotope([ABcard.center,[ABcard.generators,W.generators,V.generators,AV.generators]]);%AB * card_zono + W_sdp;
       
        
        %same as before convert R to interval
        %extract center
        c = R{i+1}.Z(:,1);
        
        %determine left and right limit
        delta = sum(abs(R{i+1}.Z),2) - abs(c);
        leftLimit{i} = c - delta;
        rightLimit{i} = c + delta;
        
        
        Constraints = [Constraints,...
            u_model{i} == U.center + U.generators * alpha_u(:,i),...
            y_model{i+1} - sinf{i} == leftLimit{i},...
            y_model{i+1} + ssup{i} == rightLimit{i},...
            y_model{i+1}   - sinf{i} >= intc.inf,...
            y_model{i+1}   + ssup{i} <= intc.sup,...
            sinf{i} >= zeros(dim_x,1),...
            ssup{i} >= zeros(dim_x,1),...
            alpha_u(i) <= ones(dim_u, 1) , ...
            alpha_u(i) >= -1*ones(dim_u, 1), ...
            ];
    end
    
    
    
    Cost_model=0;
    for i=1:N
        Cost_model = Cost_model + (y_model{i+1}-ref)'*Qy*(y_model{i+1}-ref)+ (u_model{i}-uref)'*Qu*(u_model{i}-uref);
    end
    options = sdpsettings('verbose',0,'solver','mosek');
    tic
    Problem = optimize(Constraints,Cost_model,options);
    execTimeRMPC=[execTimeRMPC,toc];
    Objective = double(Cost_model);
    uPred_model(:,timesteps) = double(u_model{1});
    YPred_model(:,timesteps+1) = double(y_model{2});

    
    % apply the optimal control input to the plant 
    w_point = randPoint(W);
    v_point = randPoint(V);
    y_t(:,timesteps+1) = sys_d.A * y_t(:,timesteps) + sys_d.B * uPred(:, timesteps) + w_point +v_point - sys_d.A *v_point;
    y_t_model(:,timesteps+1) = sys_d.A * y_t_model(:,timesteps) + sys_d.B * uPred_model(:, timesteps) + w_point +v_point - sys_d.A *v_point;
    
    yt2ref(timesteps)= norm(y_t(:,timesteps)-ref,2);
    yt2ref_model(timesteps)= norm(y_t_model(:,timesteps)-ref,2);
    halt = 1;
end



Cost_model=0;
for i=1:timesteps
    Cost_model_vec(i) = (y_t_model(:,i+1)-ref)'*Qy*(y_t_model(:,i+1)-ref)+ (uPred_model(:,i)-uref)'*Qu*(uPred_model(:,i)-uref);
    Cost_model = Cost_model + Cost_model_vec(i);

end

Cost=0;
for i=1:timesteps
    Cost_vec(i) = (y_t(:,i+1)-ref)'*Qy*(y_t(:,i+1)-ref)+ (uPred(:,i)-uref)'*Qu*(uPred(:,i)-uref);
    Cost = Cost + Cost_vec(i);
end
meanZPCtime= mean(execTimeZPC)
stdZPCtime= std(execTimeZPC)
meanRMPCtime= mean(execTimeRMPC)
stdRMPCtime= std(execTimeRMPC)

%save the workspace
save('workspaces\ZPC_cart_unstable_pi4');
%next run plotPolyZono for plotting