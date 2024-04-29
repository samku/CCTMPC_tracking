clc
clear all
addpath(genpath('/Work_stuff/IMT_set_work/softwares/yalmip_latest'))
addpath(genpath('/Work_stuff/IMT_set_work/Work_set_based/utilities'))
addpath(genpath('/Work_stuff/IMT_set_work/softwares/mpt_files'))
addpath(genpath('/Work_stuff/IMT_set_work/softwares/mosek'))
rmpath(genpath('/Work_stuff/IMT_set_work/softwares/cplex'))

nx=2;
nu=1;
A=[1.1 1; 0 1];
B=[0.5;1];
Bw=[1 0; 0 1];

HX=[eye(nx);-eye(nx)];
hX=[5;3;5;2];
mX=size(HX,1);

HU=[eye(nu);-eye(nu)];
hU=[2;1];
mU=size(HU,1);

HW=[eye(nx);-eye(nx)];
hW=[0;0.5;0;0.5];
mW=size(HW,1);

ny=1;
C=[1 0];
D=0;


%%
%Config constrained polytope for RCI set
addpath(genpath('/Work_stuff/IMT_set_work/softwares/gurobi10.0.3_linux64'))
m=16;
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
Qx_center = 10*eye(nx);
Qy_center = [V_bar -Mx]'*Qx_center*[V_bar -Mx];
Qy_extended = blkdiag(Qy_dist,zeros(n_theta,n_theta))+Qy_center;

Qu_center = eye(nu);
Qub_center = [U_bar -Mu]'*Qu_center*[U_bar -Mu];
Qub_extended =  blkdiag(Qub_dist,zeros(n_theta,n_theta))+Qub_center;

Qz_track = eye(ny);

%%
%Limon mRPI set
y=sdpvar(m,1);
us = sdpvar(m_bar*nu,1);
con=[E*y<=0];
us_r = reshape(us,nu,m_bar);
for k=1:m_bar
    con=[con; Y*(A*V{k}*y+B*us_r(:,k))+d<=y];
    con=[con; HU*us_r(:,k)<=hU; HX*V{k}*y<=hX];
end
cost = [y;zeros(nu,1)]'*Qy_extended*[y;zeros(nu,1)]+[us;zeros(nu,1)]'*Qub_extended*[us;zeros(nu,1)];
optimize(con,cost,sdpsettings('solver','gurobi'));

y_mRPI=value(y);
u_mRPI=value(us);

U_mRPI_set=minHRep(Polyhedron(value(us_r)'));
HU_mRPI = U_mRPI_set.H(:,1:end-1);
hU_mRPI = U_mRPI_set.H(:,end);

%%
%Rigid tube
Q_limon = Qx_center;
R_limon = Qu_center;
[K_limon, P_limon]=dlqr(A,B,Q_limon,R_limon);
K_limon = -K_limon;

A_aug = [A+B*K_limon B*Mu-B*K_limon*Mx; zeros(nu,nx) eye(nu)];
hX_tight = [];
for i=1:mX
    xx = cplexlp(-HX(i,:)',Y,y_mRPI);
    hX_tight(i,1)=hX(i,1)-HX(i,:)*xx;
end
hU_tight = [];
for i=1:mU
    uu = cplexlp(-HU(i,:)',HU_mRPI,hU_mRPI);
    hU_tight(i,1)=hU(i,1)-HU(i,:)*uu;
end

lambda = 0.999;
H_aug = [HX zeros(mX,nu); HU*K_limon HU*Mu-HU*K_limon*Mx; zeros(mX,nx) HX*Mx; zeros(mU,nx) HU*Mu];
h_aug = [hX_tight; hU_tight; lambda*hX_tight; lambda*hU_tight];

H_RPI = H_aug;
h_RPI = h_aug;
t=0;
RPI_found = 0;
while RPI_found==0
    t=t+1
    H_LHS = [H_RPI*A_aug; H_aug];
    h_RHS = [h_RPI; h_aug];
    step_set = minHRep(Polyhedron(H_LHS,h_RHS));
    H_RPI = step_set.H(:,1:end-1);
    h_RPI = step_set.H(:,end);
    %Verify invariance
    num_rows=size(H_RPI,1);
    lhs_vec = [];
    for i=1:num_rows
        xx = cplexlp(-(H_RPI(i,:)*A_aug)',H_RPI,h_RPI);
        lhs_vec(i,1)=H_RPI(i,:)*A_aug*xx;
    end
    diff_vec = lhs_vec-h_RPI;
    if max(diff_vec)<=1e-6
        RPI_found=1;
    end
end

%%
%Elastic tube
Q_elastic = Qx_center;
R_elastic = Qu_center;
[K_elastic, P_elastic]=dlqr(A,B,Q_elastic,R_elastic);
K_elastic = -K_elastic;

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

%Dual problems to get multipliers
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


%Terminal set with closed loop dynamics usual as limon for nominal
% and a^+ = L_RPI*a+d for RPI set
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
    %Verify invariance
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
%Compare stabilizable regions of Limon and ours
N=5;
gamma = 0.95;

c_vec = sdpvar(nx,1);

theta = sdpvar(nu,1);
x0_limon = sdpvar(nx,1);
x_limon = sdpvar(nx,N,'full');
u_limon = sdpvar(nu,N-1,'full');
con=[Y*x0_limon<=y_mRPI+Y*x_limon(:,1); HX*Mx*theta<=hX; HU*Mu*theta<=hU];
cost = -c_vec'*x0_limon;
for t=1:N-1
    con=[con; x_limon(:,t+1)==A*x_limon(:,t)+B*u_limon(:,t)];
    con=[con; HX*x_limon(:,t)<=hX_tight];
    con=[con; HU*u_limon(:,t)<=hU_tight];
end
con=[con; H_RPI*[x_limon(:,N);theta]<=h_RPI];

theta_elastic = sdpvar(nu,1);
x0_elastic = sdpvar(nx,1);
x_elastic = sdpvar(nx,N,'full');
u_elastic = sdpvar(nu,N-1,'full');
a_elastic = sdpvar(m,N,'full');
con=[con; Y*(x0_elastic-x_elastic(:,1))<=a_elastic(:,1); HX*Mx*theta_elastic<=hX; HU*Mu*theta_elastic<=hU];
cost = cost-c_vec'*x0_elastic;
for t=1:N-1
    con=[con; Y*(A*x_elastic(:,t)+B*u_elastic(:,t))+L_RPI*a_elastic(:,t)+d<=Y*x_elastic(:,t+1)+a_elastic(:,t+1)];
    con=[con; HX*x_elastic(:,t)+L_HX*a_elastic(:,t)<=hX];
    con=[con; HU*u_elastic(:,t)+L_HU*a_elastic(:,t)<=hU];
    con=[con; a_elastic(:,t)>=0];
end
con=[con; H_RPI_elastic*[x_elastic(:,N);theta;a_elastic(:,N)]<=h_RPI_elastic];

y=sdpvar(m,N,'full');
u=sdpvar(nu*m_bar,N,'full');
x0_CCTMPC=sdpvar(nx,1);
ys=sdpvar(m,1);
us=sdpvar(m_bar*nu,1);
z_des=sdpvar(ny,1);
con=[con; Y*x0_CCTMPC<=y(:,1); E*y<=repmat(zeros(size(E,1),1),1,N); E*ys<=zeros(size(E,1),1)];
us_r = reshape(us,nu,m_bar);
for k=1:m_bar
    con=[con; Y*(A*V{k}*ys+B*us_r(:,k))+d<=ys];
    con=[con; HU*us_r(:,k)<=hU; HX*V{k}*ys<=hX];
end
cost = cost-c_vec'*x0_CCTMPC;
for t=1:N-1
    u_loc=reshape(u(:,t),nu,m_bar);
    for k=1:m_bar
        con=[con; HX*V{k}*y(:,t)<=hX; HU*u_loc(:,k)<=hU];
        con=[con; Y*(A*V{k}*y(:,t)+B*u_loc(:,k))+d<=y(:,t+1)];
    end
end
u_loc=reshape(u(:,N),nu,m_bar);
for k=1:m_bar
    con=[con; HX*V{k}*y(:,N)<=hX; HU*u_loc(:,k)<=hU];
    con=[con; Y*(A*V{k}*y(:,N)+B*u_loc(:,k))+d<=gamma*y(:,N)+(1-gamma)*ys];
end

get_stab_reg = optimizer(con,cost,sdpsettings('solver','gurobi'),{c_vec},{x0_limon,x0_elastic,x0_CCTMPC});

H_stab=get_param_matrices_all_permutations(nx,nx,200,200);
h_limon = [];
x0_limon = [];
h_elastic = [];
x0_elastic = [];
h_CCTMPC = [];
x0_CCTMPC = [];
for i=1:size(H_stab,1)
    i
    soln = get_stab_reg{H_stab(i,:)'};
    x0_limon(:,i)=soln{1};
    x0_elastic(:,i)=soln{2};
    x0_CCTMPC(:,i)=soln{3};
    h_limon(i,1)=H_stab(i,:)*x0_limon(:,i);
    h_elastic(i,1)=H_stab(i,:)*x0_elastic(:,i);
    h_CCTMPC(i,1)=H_stab(i,:)*x0_CCTMPC(:,i);
end

figure
polyopts.alpha=0.4;
polyopts.color='red';
plot(polytope(H_stab,h_CCTMPC),polyopts)
hold on
polyopts.alpha=0.4;
polyopts.color='green';
plot(polytope(H_stab,h_elastic),polyopts)
hold on
polyopts.alpha=0.4;
polyopts.color='blue';
plot(polytope(H_stab,h_limon),polyopts)

polyopts.color='green';
polyopts.alpha=0.0001;
H_MRCI = HX;
h_MRCI = hX;
for i=1:100
    i
    d_loc = [];
    for j=1:size(H_MRCI,1)
        ww = cplexlp(-(H_MRCI(j,:)*Bw)',HW,hW);
        d_loc(j,1)=H_MRCI(j,:)*Bw*ww;
    end
    H_LHS = [H_MRCI*A H_MRCI*B; HX zeros(mX,nu); zeros(mU,nx) HU];
    h_RHS = [h_MRCI-d_loc; hX; hU];
    set_full = minHRep(projection(Polyhedron(H_LHS,h_RHS),1:nx));
    H_MRCI = set_full.H(:,1:end-1);
    h_MRCI = set_full.H(:,end);
    plot(polytope(H_MRCI,h_MRCI),polyopts)
    hold on
    pause(0.1)
    
end

%%
vert_MRCI=Polyhedron(H_MRCI,h_MRCI).V';
x_pts = sdpvar(nx,size(vert_MRCI,2));
b_pts = sdpvar(nx,size(vert_MRCI,2));
epsilon = sdpvar(1,1);
con=[vert_MRCI==x_pts+b_pts; H_stab*x_pts<=repmat(h_CCTMPC,1,size(vert_MRCI,2)); [eye(nx);-eye(nx)]*b_pts<=repmat(epsilon*ones(2*nx,1),1,size(vert_MRCI,2))];
optimize(con,epsilon);
epsilon_CCTMPC=value(epsilon)

x_pts = sdpvar(nx,size(vert_MRCI,2));
b_pts = sdpvar(nx,size(vert_MRCI,2));
epsilon = sdpvar(1,1);
con=[vert_MRCI==x_pts+b_pts; H_stab*x_pts<=repmat(h_limon,1,size(vert_MRCI,2)); [eye(nx);-eye(nx)]*b_pts<=repmat(epsilon*ones(2*nx,1),1,size(vert_MRCI,2))];
optimize(con,epsilon);
epsilon_limon=value(epsilon)

x_pts = sdpvar(nx,size(vert_MRCI,2));
b_pts = sdpvar(nx,size(vert_MRCI,2));
epsilon = sdpvar(1,1);
con=[vert_MRCI==x_pts+b_pts; H_stab*x_pts<=repmat(h_elastic,1,size(vert_MRCI,2)); [eye(nx);-eye(nx)]*b_pts<=repmat(epsilon*ones(2*nx,1),1,size(vert_MRCI,2))];
optimize(con,epsilon);
epsilon_elastic=value(epsilon)