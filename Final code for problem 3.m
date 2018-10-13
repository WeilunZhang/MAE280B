%please make sure to run by section to obtain answers, thanks!
%%
%problem3,a)
clc
clear all

%define state space
A = [0 0 1 0;0 0 0 1;-2 1 0 0;1 -1 0 0];
Bu = [0 0;0 0;1 0;0 1];
Bw =[0 0;0 0;0 0;0 1];
n = size(A,1);
m = size(Bu,2);
r = size(Bw,2);
Cz=[0 0 0 0; 0 0 0 0;1 0 0 0; 0 1 0 0];
Dzu=[1 0; 0 1; 0 0; 0 0];
Cz0 = Cz(1:2,:);
Cz1 = Cz(3,:);
Cz2 = Cz(4,:);
Dzu0 = Dzu(1:2,:);
Dzu1 = Dzu(3,:);
Dzu2 = Dzu(4,:);
miu=0.5;
W =eye(r);

%define LMI variables
tol=1e-6;
X = sdpvar(n,n,'symmetric');
Z0 = sdpvar(2,2);
Z1 = sdpvar(1,1);
Z2 = sdpvar(1,1);
L = sdpvar(m,n);

%define LMI
LMI1 = A*X+X*A'+Bu*L+L'*Bu'+Bw*W*Bw';
LMI2 = [Z0, (Cz0*X+Dzu0*L); (X*Cz0'+L'*Dzu0'), X];
LMI3 = X;
LMI4 = [Z1, (Cz1*X+Dzu1*L); (X*Cz1'+L'*Dzu1'), X];
LMI5 = [Z2, (Cz2*X+Dzu2*L); (X*Cz2'+L'*Dzu2'), X];
H=[1 0;-1 0;0 1;0 -1];
limit=[tol;-miu^2;tol;-miu^2];
LMI6 = H*[trace(Z1);trace(Z2)]-limit;
LMI = [LMI1 <= -tol,LMI2 >= tol,LMI3>=tol,LMI4>=tol,LMI5>=tol,LMI6>=tol];
options = sdpsettings('solver','sedumi');
solution =solvesdp(LMI,trace(Z0),options);

%============================================================================
K = double(L)/double(X);
Z0=trace(double(Z0));
J0=trace((Cz0*double(X)+Dzu0*double(L))*inv(double(X))*(Cz0*double(X)+Dzu0*double(L))')
Z1=double(Z1)
J1=(Cz1*double(X)+Dzu1*double(L))*inv(double(X))*(Cz1*double(X)+Dzu1*double(L))'
Z2=double(Z2)
J2=(Cz2*double(X)+Dzu2*double(L))*inv(double(X))*(Cz2*double(X)+Dzu2*double(L))'

%define cloop
Acl = A + Bu*K
Ccl = Cz0 + Dzu0*K;
%eig(Acl)
syscl=ss(Acl,Bu,Ccl,Dzu0);
step(syscl)
%obtain controller
sysctrl=ss(A,Bu,K,0);
tfctrl=tf(sysctrl)

%%
%problem3,d)m1&wall
clc
clear all
%define system
A = [0 1;-2 0];
Bu = [0 0;1 1];
Bw =[0 0;0 0];
n = size(A,1);
m = size(Bu,2);
r = size(Bw,2);
Cz0 = [0 0];
Cz1 = [1 0];
Dzu0 = [1 0];
Dzu1 = [0 0];
miu=0.5;
W =eye(r);

%define LMI variables
tol=1e-5;
X = sdpvar(n,n,'symmetric');
Z0 = sdpvar(1,1);
Z1 = sdpvar(1,1);
L = sdpvar(m,n);

%formulate LMI
LMI1 = A*X+X*A'+Bu*L+L'*Bu'+Bw*W*Bw';
LMI2 = [Z0, (Cz0*X+Dzu0*L); (X*Cz0'+L'*Dzu0'), X];%J0
LMI3 = X;
LMI4 = [Z1, (Cz1*X+Dzu1*L); (X*Cz1'+L'*Dzu1'), X];%J1
TM=[1 0;-1 0];
lim=[tol;-miu^2];
LMI6 = TM*[trace(Z1);0]-lim;
LMI = [LMI1 <= -tol,LMI2 >= tol,LMI3>=tol,LMI4>=tol,LMI6>=tol];
options = sdpsettings('solver','sedumi');
solution =solvesdp(LMI,trace(Z0),options);

K = double(L)/double(X);
Acl = A + Bu*K
%eig(Acl)
%%
%problem3,d)m1&m2
clc
clear all
%define system
A = [0 1;-1 0];
Bu = [0 0;1 1];
Bw =[0 0;0 1];
n = size(A,1);
m = size(Bu,2);
r = size(Bw,2);
Cz0 = [0 0];
Cz1 = [1 0];
Dzu0 = [1 0];
Dzu1 = [0 0];
miu=0.5;
W =eye(r);

%define LMI variables
%tol=2.6e-3;
tol=1e-5;
X = sdpvar(n,n,'symmetric');
Z0 = sdpvar(1,1);
Z1 = sdpvar(1,1);
L = sdpvar(m,n);
%formulate LMI
LMI1 = A*X+X*A'+Bu*L+L'*Bu'+Bw*W*Bw';
LMI2 = [Z0, (Cz0*X+Dzu0*L); (X*Cz0'+L'*Dzu0'), X];%J0
LMI3 = X;
LMI4 = [Z1, (Cz1*X+Dzu1*L); (X*Cz1'+L'*Dzu1'), X];%J1
TM=[1 0;-1 0];
lim=[tol;-miu^2];
LMI6 = TM*[trace(Z1);0]-lim;
LMI = [LMI1 <= -tol,LMI2 >= tol,LMI3>=tol,LMI4>=tol,LMI6>=tol,abs(trace(Z1)-miu^2)<=tol];
options = sdpsettings('solver','sedumi');
solution =solvesdp(LMI,trace(Z0),options);


K = double(L)/double(X)
Acl = A + Bu*K
%eig(Acl)





