clc
clear all
addpath(genpath('/Work_stuff/IMT_set_work/softwares/yalmip_latest'))
addpath(genpath('/Work_stuff/IMT_set_work/Work_set_based/utilities'))
addpath(genpath('/Work_stuff/IMT_set_work/softwares/mpt_files'))
addpath(genpath('/Work_stuff/IMT_set_work/softwares/mosek'))
rmpath(genpath('/Work_stuff/IMT_set_work/softwares/cplex'))


%%
nx=2;
nu=1;
A=[1.1 1; 0 1];
B=[0.5;1];
Bw=[1 0; 0 1];

HX=[eye(nx);-eye(nx)];
hX=[5;3; 5;2];
mX=size(HX,1);

HU=[eye(nu);-eye(nu)];
hU=[2; 1];
mU=size(HU,1);

HW=[eye(nx);-eye(nx)];
hW=[0;0.5; 0;0.5];
mW=size(HW,1);

ny=1;
C=[1 0];
D=0;

%%
%Config constrained polytope for RCI set
addpath(genpath('/Work_stuff/IMT_set_work/softwares/gurobi10.0.3_linux64'))
m=22;
Y=get_param_matrices_all_permutations(nx,nx,m,m);
m=size(Y,1);

y0=ones(m,1);
WH_P=Polyhedron(Y,y0);
WH_vert=WH_P.V';

%Get index sets
WC=nchoosek(1:size(Y,1),nx);
V={};
index_mat={};
one_Vi={};
for i=1:size(WH_vert,2)
    for j=1:size(WC,1)
        V_test=Y(WC(j,:),:);
        h_test=y0(WC(j,:));
        if norm(V_test*WH_vert(:,i)-h_test,'inf')<=1e-9
            index_mat{i}=WC(j,:)';
            ones_mat=zeros(nx,size(Y,1));
            for k=1:nx
                ones_mat(k,WC(j,k))=1;
            end
            one_Vi{i}=ones_mat;
            V{i}=inv(V_test)*one_Vi{i};
        end
    end
end

E=[];
for i=1:size(WH_vert,2)
    E=[E; Y*V{i}-eye(m)];
end
E=sparse(E);
m_bar=length(V); %Number of vertices
l=size(E,1);

addpath(genpath('/Work_stuff/IMT_set_work/softwares/cplex'))
d=[];
for i=1:m
    ww=cplexlp(-(Y(i,:)*Bw)',HW,hW);
    d(i,1)=Y(i,:)*Bw*ww;
end

%%
Qx_dist=10*eye(nx);
Qu_dist=eye(nu);

V_bar=zeros(size(V{1},1),size(V{1},2));
for i=1:length(V)
    V_bar=V_bar+V{i};
end
V_bar=V_bar/length(V);

Qy_dist=zeros(m,m);
for i=1:length(V)
    Qy_dist=Qy_dist+(V{i}-V_bar)'*Qx_dist*(V{i}-V_bar);
end
Iu_mats={};
U_bar=0;
for i=1:length(V)
    Iu_mats{i}=sparse(nu,nu*length(V));
    Iu_mats{i}(:,(i-1)*nu+1:i*nu)=eye(nu);
    U_bar=U_bar+Iu_mats{i};
end
Qub_dist=zeros(nu*m_bar,nu*m_bar);
U_bar=U_bar/m_bar;
for i=1:m_bar
    Qub_dist=Qub_dist+(Iu_mats{i}-U_bar)'*Qu_dist*(Iu_mats{i}-U_bar);
end

M = null([A-eye(nx) B]);
n_theta = size(M,2);
Mx = M(1:nx,:);
Mu = M(nx+1:end,:);
Qx_center = eye(nx);
Qy_center = [V_bar -Mx]'*Qx_center*[V_bar -Mx];
Qy_extended = blkdiag(Qy_dist,zeros(n_theta,n_theta))+Qy_center;

Qu_center = eye(nu);
Qub_center = [U_bar -Mu]'*Qu_center*[U_bar -Mu];
Qub_extended =  blkdiag(Qub_dist,zeros(n_theta,n_theta))+Qub_center;

Qz_track = eye(ny);
Qtheta_track = (C*Mx)'*Qz_track*(C*Mx);
Qy_extended=Qy_extended+blkdiag(zeros(m,m),Qtheta_track);

%%
Q_elastic = Qx_center;
R_elastic = Qu_center;
[K_elastic, P_elastic]=dlqr(A,B,Q_elastic,R_elastic);
K_elastic = -K_elastic;


%%
%Rakovic mRPI set
y=sdpvar(m,1);
us = sdpvar(nu,m,'full');
xs = sdpvar(nx,m,'full');
c_vec = sdpvar(m,1);
con=[];
for i=1:m
    con=[con; c_vec(i)<=Y(i,:)*(A+B*K_elastic)*xs(:,i)];
    con=[con; Y*xs(:,i)<=c_vec+d];
end
cost = -sum(c_vec);
optimize(con,cost,sdpsettings('solver','gurobi'));

c_vec = value(c_vec);
y_mRPI_elastic=c_vec+d;

c_ver = [];
for i=1:m
    xx = cplexlp(-(Y(i,:)*(A+B*K_elastic))',Y,y_mRPI_elastic);
    c_ver(i,1)=Y(i,:)*(A+B*K_elastic)*xx;
end

L_RPI = [];
L_HX = [];
L_HU = [];
for i=1:m
    lx = cplexlp(y_mRPI_elastic,-eye(m),zeros(m,1),Y',(Y(i,:)*(A+B*K_elastic))');
    L_RPI(i,:)=lx';
end
for i=1:mX
    lx = cplexlp(y_mRPI_elastic,-eye(m),zeros(m,1),Y',HX(i,:)');
    L_HX(i,:)=lx';
end
for i=1:mU
    lx = cplexlp(y_mRPI_elastic,-eye(m),zeros(m,1),Y',K_elastic'*HU(i,:)');
    L_HU(i,:)=lx';
end    
   
Q_RPI = eye(m);
P_RPI = sdpvar(m,m);
LMI_val = sdpvar(m,m);
con=[LMI_val==-Q_RPI+P_RPI-L_RPI'*P_RPI*L_RPI; LMI_val>=0; P_RPI>=0];
optimize(con,trace(P_RPI),sdpsettings('solver','mosek'));
P_RPI=value(P_RPI);

%%
%Terminal set - dynamics of [z;theta;a]
A_aug_elastic = [A+B*K_elastic B*Mu-B*K_elastic*Mx zeros(nx,m); ...
             zeros(nu,nx)    eye(nu)    zeros(nu,m); ...
             zeros(m,nx)     zeros(m,nu)  L_RPI];
d_aug_elastic = [zeros(nx+nu,1);d];


lambda = 0.999;
H_aug_elastic = [HX zeros(mX,nu) L_HX; ...
         HU*K_elastic HU*Mu-HU*K_elastic*Mx L_HU; ...
         zeros(mX,nx) HX*Mx  zeros(mX,m); ... 
         zeros(mU,nx) HU*Mu  zeros(mU,m); ...
         zeros(m,nx)  zeros(m,nu) -eye(m)];
     
h_aug_elastic = [hX; hU; lambda*hX; lambda*hU; zeros(m,1)];

H_RPI_elastic = H_aug_elastic;
h_RPI_elastic = h_aug_elastic;
RPI_found = 0;
t=0;
while RPI_found==0
    t=t+1
    H_RPI_elastic = [H_RPI_elastic; H_aug_elastic*A_aug_elastic^t];
    rhs_reduce = zeros(size(h_aug_elastic,1),1);
    for jj=t-1:-1:0
        rhs_reduce=rhs_reduce+H_aug_elastic*(A_aug_elastic^jj)*d_aug_elastic;
    end
    h_RPI_elastic = [h_RPI_elastic; h_aug_elastic-rhs_reduce];
    %Verify invariance H(Ax+d)<=h for all x:Hx<=h
    num_rows=size(H_RPI_elastic,1);
    lhs_vec = [];
    for i=1:num_rows
        xx = cplexlp(-(H_RPI_elastic(i,:)*A_aug_elastic)',H_RPI_elastic,h_RPI_elastic);
        lhs_vec(i,1)=H_RPI_elastic(i,:)*A_aug_elastic*xx+H_RPI_elastic(i,:)*d_aug_elastic;
    end
    diff_vec = lhs_vec-h_RPI_elastic;
    if max(diff_vec)<=1e-6
        RPI_found=1;
    end
end

%%
N=10;
z_des = sdpvar(ny,1);
theta = sdpvar(nu,1);
x0 = sdpvar(nx,1);
x = sdpvar(nx,N,'full');
u = sdpvar(nu,N-1,'full');
a = sdpvar(m,N,'full');
con=[Y*(x0-x(:,1))<=a(:,1); HX*Mx*theta<=hX; HU*Mu*theta<=hU];
cost = (z_des-[C D]*M*theta)'*Qz_track*(z_des-[C D]*M*theta);
for t=1:N-1
    con=[con; Y*(A*x(:,t)+B*u(:,t))+L_RPI*a(:,t)+d<=Y*x(:,t+1)+a(:,t+1)];
    con=[con; HX*x(:,t)+L_HX*a(:,t)<=hX];
    con=[con; HU*u(:,t)+L_HU*a(:,t)<=hU];
    con=[con; a(:,t)>=0];
    cost=cost+(x(:,t)-Mx*theta)'*Q_elastic*(x(:,t)-Mx*theta)+(u(:,t)-Mu*theta)'*R_elastic*(u(:,t)-Mu*theta);
    cost=cost+(a(:,t)-y_mRPI_elastic)'*Q_RPI*(a(:,t)-y_mRPI_elastic);
end
cost=cost+(x(:,N)-Mx*theta)'*P_elastic*(x(:,N)-Mx*theta)+(a(:,N)-y_mRPI_elastic)'*P_RPI*(a(:,N)-y_mRPI_elastic);
con=[con; H_RPI_elastic*[x(:,N);theta;a(:,N)]<=h_RPI_elastic];
MPC_elastic = optimizer(con,cost,sdpsettings('solver','gurobi'),{x0,z_des},{x,u,theta,a,cost});


theta = sdpvar(nu,1);
x0_des = sdpvar(nx,1);
x0 = sdpvar(nx,1);
x = sdpvar(nx,N,'full');
u = sdpvar(nu,N-1,'full');
a = sdpvar(m,N,'full');
con=[Y*(x0-x(:,1))<=a(:,1); HX*Mx*theta<=hX; HU*Mu*theta<=hU];
cost = (x0_des-x0)'*(x0_des-x0);
for t=1:N-1
    con=[con; Y*(A*x(:,t)+B*u(:,t))+L_RPI*a(:,t)+d<=Y*x(:,t+1)+a(:,t+1)];
    con=[con; HX*x(:,t)+L_HX*a(:,t)<=hX];
    con=[con; HU*u(:,t)+L_HU*a(:,t)<=hU];
    con=[con; a(:,t)>=0];
end
con=[con; H_RPI_elastic*[x(:,N);theta;a(:,N)]<=h_RPI_elastic];
MPC_elastic_projector = optimizer(con,cost,sdpsettings('solver','gurobi'),{x0_des},{x0});

%%
K_des = sdpvar(nu,nx,'full');
K_projected = sdpvar(nu,nx,'full');
con=[L_RPI*Y==Y*(A+B*K_projected); L_HU*Y==HU*K_projected];
cost=norm(K_projected-K_des,1);
K_computer = optimizer(con,cost,sdpsettings('solver','gurobi'),{K_des},{K_projected});

%%
figure
X_vert=Polyhedron(HX,hX).V';
x_sys=[-3.321;3];
x_sys=MPC_elastic_projector{x_sys};
u_sys=[];
w_sys=[];

x_MPC={};
u_MPC={};
thetas_MPC = [];
a_MPC = {};
cost_MPC=[];
theta_opt = [];

z_des=[];

W_vert=Polyhedron(HW,hW).V';

polyopts1.alpha=0.1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
polyopts1.color='red';

polyopts2.alpha=0.1;
polyopts2.color='green';

polyopts3.alpha=0.1;
polyopts3.color='black';

polyopts4.alpha=0.2;
polyopts4.color='blue';

polyopts5.alpha=0.2;
polyopts5.color='yellow';
subplot(2,2,1)
hold on

for t=1:39
    if t==1 || rem(t,20)==0
        if t<40
            z_des(:,t)=5
        else
            z_des(:,t)=-5;
        end
    else
        z_des(:,t)=z_des(:,t-1);
    end
    
    tic
    MPC_soln=MPC_elastic{x_sys(:,t),z_des(:,t)};
    toc
    x_MPC{t}=MPC_soln{1};
    u_MPC{t}=MPC_soln{2};
    thetas_MPC(:,t)=MPC_soln{3};
    a_MPC{t}=MPC_soln{4};
    cost_MPC(t)=MPC_soln{5};
    
    w_sys(:,t)=W_vert(:,randi(size(W_vert,2)));
    
    K_RPI = K_computer{10*randn(nu,nx)};    
    u_sys(:,t)=u_MPC{t}(:,1)+K_RPI*(x_sys(:,t)-x_MPC{t}(:,1));
    x_sys(:,t+1)=A*x_sys(:,t)+B*u_sys(:,t)+Bw*w_sys(:,t);

    y_curr=a_MPC{t}(:,1)+Y*x_MPC{t}(:,1);
    y_next=a_MPC{t}(:,2)+Y*x_MPC{t}(:,2);
    
    subplot(2,2,1)
    plot(polytope(Y,y_mRPI_elastic+Y*Mx*thetas_MPC(:,t)),polyopts4)
    hold on
    plot(polytope(Y,y_curr),polyopts1)
    hold on
    plot(polytope(Y,y_next),polyopts2)
    hold on
    %for i=3:N
    %    y_locc=a_MPC{t}(:,i)+Y*x_MPC{t}(:,i);
    %    plot(polytope(Y,y_locc),polyopts3)
    %    hold on
    %end
    plot(x_sys(1,:),x_sys(2,:),'black')
    hold on
    scatter(x_sys(1,t),x_sys(2,t),'redo')
    
   
    subplot(2,2,3)
    scatter(t*ones(ny,1),C*x_sys(:,t),'black*')
    hold on
    scatter(t*ones(ny,1),z_des(:,t),'redo')
    hold on

    
    subplot(2,2,4)
    semilogy(t,cost_MPC(t),'black*')
    hold on

    
    
    pause(0.0001);
 
end

