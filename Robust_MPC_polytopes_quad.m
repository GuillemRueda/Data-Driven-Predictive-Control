% ZPC: run RMPC with model using polytopes
% 
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

% Author:       Yvonne Stürz 
% Written:      25-March-2021 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


clear all
close all
rand('seed',4500);

% dimension of x
dim_x = 6;
dim_u = 3;

% Nominal system
m = 1;
I = 0.5;
l = 0.25;
g= 9.81;
u1_t = m*g/2;
u2_t = u1_t;

% System in cont time
A = [0, 0, 0, 1, 0, 0; 0, 0, 0, 0, 1, 0; 0, 0, 0, 0, 0, 1;
    0, 0, -(u1_t + u2_t)/m, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0];
B_ss = [0, 0, 0; 0, 0, 0; 0, 0, 0; 0, 0, 0; 1/m, 1/m, -1; l/I, -l/I, 0];
C = eye(dim_x);
D = zeros(dim_x, dim_u);
% define continuous time system
sys_c = ss(A,B_ss,C,D);
% convert to discrete system
samplingtime = 0.05;
sys_d = c2d(sys_c,samplingtime);
A = sys_d.A;
B = sys_d.B;
C = eye(dim_x);
E = eye(dim_x);

wfac=0.01;
vfac = 0.002;
Qy = 1e3*eye(dim_x); %define output cost matrix
Qu = 0.001*eye(dim_u);%control cost matrix

uref = [u1_t; u2_t; g];  
% ref = inv(eye(5)-sys_d.A)*sys_d.B*uref;
ref = [0; 0; 0; 0; 0; 0];

y_lb = [-10;2;-0.1;-10;-10;-10]; 
y_ub = [10;10;0.1;10;10;10]; 
intc = interval(y_lb,y_ub);


% Open loop mini-max solution

N = 2;
U = sdpvar(N,dim_u);
W = sdpvar(N,dim_x);
V = sdpvar(N,dim_x);
x = sdpvar(dim_x,1);

Vw = wfac*ones(dim_x,1);
Vw = [Vw'; -Vw'];
Pw = Polyhedron('V', Vw);
%plot(P1)
Pw.minHRep;
Pw.H;
Pw.He;

test = [0;0;0;0;0;0];
W_Constraints = [Pw.A * test <= Pw.b; Pw.Ae * test == Pw.be];


Vv = vfac*ones(dim_x, 1);
Vv = [Vv'; -Vv'];
Pv = Polyhedron('V', Vv);
%plot(P1)
Pv.minHRep;
Pv.H;
Pv.He;

test = [0;0;0;0;0;0];
V_Constraints = [Pv.A * test <= Pv.b; Pv.Ae * test == Pv.be];


Y = [];
xk = x;
for k = 1:N
 xk = A*xk + B*U(k,:)'+E*W(k,:)';
 Y = [Y; C*xk + V(k,:)'];
end

% TODO: fix this 
F = [kron(ones(N,1),y_lb) <= Y <= kron(ones(N,1),y_ub), kron(ones(N,dim_u),-32) <= U <= kron(ones(N,dim_u),46),
    kron(ones(N,1),-) <= U <= kron(ones(N,1),46),
    kron(ones(N,1),9.81) <= U <= kron(ones(N,1),9.81)];
objective = norm(Y-kron(ones(N,1),ref),2)*Qy(1) + norm(U-kron(ones(N,dim_u),uref),2)*Qu(1);

G = []; 
for k = 1:N
    G = [G, blkdiag(Pw.A,Pv.A) * [W(k,:)'; V(k,:)'] <= [Pw.b; Pv.b], blkdiag(Pw.Ae,Pv.Ae) * [W(k,:)'; V(k,:)'] == [Pw.be;Pv.be]];
end


[Frobust,h] = robustify(F + G,objective,[],[W;V]);


xk = y0;
uk = [];
Y = y0;
ops = sdpsettings;
maxsteps = 80;
execTime=[];
for i = 1:maxsteps
    tic
    optimize([Frobust, x == xk(:,end)],h,ops);
    execTime=[execTime toc];
    xk = [xk A*xk(:,end) + B*value(U(1)) + E*wfac*(-1+2*rand(1)*ones(5,1))];
    Y = [Y, C*xk(:,end) + vfac*(-1+2*rand(1)*ones(5,1))];
    uk = [uk value(U(1))];
end

Cost_rob_ol_tot=0;
Cost_rob_ol=[];
for i = 1:maxsteps
    Cost_rob_ol = [Cost_rob_ol, (Y(:,i+1)-ref)'*Qy*(Y(:,i+1)-ref)+ (uk(:,i)-uref)'*Qu*(uk(:,i)-uref)];
    Cost_rob_ol_tot = Cost_rob_ol_tot + (Y(:,i+1)-ref)'*Qy*(Y(:,i+1)-ref)+ (uk(:,i)-uref)'*Qu*(uk(:,i)-uref);
    yt2ref_poly(i) = norm(Y(:,i)-ref,2);
end
Cost_rob_ol_tot
meanExecTime=mean(execTime)
stdExecTime= std(execTime)
% figure(1)
% %plot([C*xk  + vfac*(-1+2*rand(1)*ones(5,1))]')
% plot([Y]')
% hold on, plot(kron(ones(1,100),ref)')
% figure(2)
% %plot([C*xk  + vfac*(-1+2*rand(1)*ones(5,1))]')
% plot([Y]')
% hold on, plot(kron(ones(1,100),ref)')
% hold on, plot(kron(ones(1,100),y_lb)') 
% hold on, plot(kron(ones(1,100),y_ub)') 

save('workspaces\poly');

