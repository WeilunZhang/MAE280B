%Please make sure to run by section to get the answers, thanks!
%%
%problem1a)1.2.3
clc 
clear all
%define the system
%==================================================================
A=[-.0558 -.9968 .0802 .0415;
     .598 -.115 -.0318 0;
    -3.05 .388 -.4650 0;
        0 0.0805 1 0];
Bu=[ .00729  0;
   -0.475   0.00775;
    0.153   0.143;
     0      0];
Bw=[0.0075 0.0026;0.0081 0.0066;0.0050 0.0075;0.0040 0.0046];
%rank([A,Bu]);
Cy=[0 1 0 0];%measure x2
%Cy=[0 0 0 1];%measure x4
%rank([A;Cy]);
Dyw=[0.0067 0.0089];
Dyu=Cy*Bu(:,1);
Cz=[eye(4);zeros(1,4)];
Dzu=[0 0 0 0 1]';
W=eye(2);
%1.a)2.=========================================================

Q1=Cz'*Cz;
Q2=Bw*W*Bw';
R1=Dzu'*Dzu;
R2=Dyw*W*Dyw';
[X,L1,G1] = care(A,Bu(:,1),Q1,R1);
[Y,L2,G2] = care(A',Cy',Q2,R2);
K=-inv(Dzu'*Dzu)*Bu(:,1)'*X;
K=[K;zeros(1,4)];
F=-Y*Cy'*inv(Dyw*W*Dyw');
F=[F,zeros(4,1)];

%1,a)3.==========================================================
%controller transfer function and state space
Cyt=[0 1 0 0; 0 0 0 1];
Actrl=A+Bu*K+F*Cyt;
Bctrl=-F;
Cctrl=K;
Dctrl=0;
sysctrl=ss(Actrl,Bctrl,Cctrl,Dctrl)
tfctrl=tf(sysctrl)

%%
%problem1a)4 .5

A=[-.0558 -.9968 .0802 .0415;
     .598 -.115 -.0318 0;
    -3.05 .388 -.4650 0;
        0 0.0805 1 0];
Bu=[ .00729  0;
   -0.475   0.00775;
    0.153   0.143;
     0      0];
Bw=[0.0075 0.0026;0.0081 0.0066;0.0050 0.0075;0.0040 0.0046];
Dyw=[0.0067 0.0089];
 %Bw =[0.0042 0.0066;0.0092 0.0004;0.0079 0.0085; 0.0096 0.0093];
Cy=[0 1 0 0];%measure x2
 %Cy=[0 0 0 1];%measure x4
% Dyw=[0.0068 0.0076];
Dyu=Cy*Bu(:,1);
Cz=[eye(4);zeros(1,4)];
Dzu=[0 0 0 0 1]';
W=eye(2);
Q1=Cz'*Cz;
Q2=Bw*W*Bw';
R1=Dzu'*Dzu;
R2=Dyw*W*Dyw';
[X,L1,G1] = care(A,Bu(:,1),Q1,R1);
[Y,L2,G2] = care(A',Cy',Q2,R2);
K=-inv(Dzu'*Dzu)*Bu(:,1)'*X;
K=[K;zeros(1,4)];
F=-Y*Cy'*inv(Dyw*W*Dyw');
F=[F,zeros(4,1)];

Cyt=[0 1 0 0; 0 0 0 1];
Actrl=A+Bu*K+F*Cyt;
Bctrl=-F;
Cctrl=K;
Dctrl=0;
sysctrl=ss(Actrl,Bctrl,Cctrl,Dctrl)
tfctrl=tf(sysctrl)
Aoloop=A;
Boloop=Bu;
Coloop=Cyt;
Doloop=0;
sysoloop=ss(Aoloop,Boloop,Coloop,Doloop);

% Acloop=[A Bu(:,1)*K;-F*Cyt A+Bu(:,1)*K+F*Cyt];
% Bcloop=[Bw;-F*Dyw];
% Ccloop=[Cz Dzu*K];
% Dcloop=0;
syscloop=feedback(sysoloop,sysctrl,1)
% syscloop=ss(Acloop,Bcloop,Ccloop,Dcloop);
axis(gca,'normal')
h1 = pzplot(syscloop);
[Wn,zeta] = damp(syscloop);%coefficients satisfied!

%problem1a)5
impulseplot(sysoloop,'b--',syscloop,'r',20)
%figure(1,1)corresponds to 'rudder to yaw rate'
%figure(2,2)corresponds to 'aileron to bank angle'
%if we want to use this same cheap controller to control also aileron to bank
%angle we only need to change Cy to Cyt

%We can also use the same method to compute a cheap controller based on the second column of Bu...
%in that case for aileron to bank angle just change Bu(:,1) to Bu(1,:) 

%please note that we can use this some procedures to design a controller
%based on both u1 and u2, however, that controller will not directly show
%how u1 change x2 and u2 change x4...
%%
clc 
clear all
%problem1c)
%reconstruct the state space
A=[-.0558 -.9968 .0802 .0415;
     .598 -.115 -.0318 0;
    -3.05 .388 -.4650 0;
        0 0.0805 1 0];
Bu=[ .00729  0;
   -0.475   0.00775;
    0.153   0.143;
     0      0];
Buk=A*Bu;
Bw=[0.0066 0.0033;0.0025 0.0089;0.0076 0.0046;0.0082 0.0038];
Dyw=[0.0087 0.0057];
Cy=[0 1 0 0];
Cyt=[0 1 0 0;0 0 0 1];
%rank([A;Cy]);
Dyu=Cyt*Bu;
Cz=[eye(4);zeros(1,4)];
Dzu=[0 0 0 0 1]';
W=eye(2);

Aoloop=A;
Boloop=Bu;
Coloop=Cyt;
Doloop=0;
sysoloop=ss(Aoloop,Boloop,Coloop,Doloop);

Q1=Cz'*Cz;
Q2=Bw*W*Bw';
R1=Dzu'*Dzu;
R2=Dyw*W*Dyw';
[X,L1,G1] = care(A,Buk(:,1),Q1,R1);
[Y,L2,G2] = care(A',Cy',Q2,R2);
K=-inv(Dzu'*Dzu)*Buk(:,1)'*X;
F=-Y*Cy'*inv(Dyw*W*Dyw');
K=[K;zeros(1,4)];
F=[F,zeros(4,1)];
Actrl=A+Buk*K+F*Cyt;
Bctrl=-F;
Cctrl=K;
Dctrl=0;
sysctrl=ss(Actrl,Bctrl,Cctrl,Dctrl);
tfctrl=tf(sysctrl);%C(s)tilta
Cs=ss(Actrl-Bctrl*Dyu*Cctrl,Bctrl,Cctrl,Dctrl);
tfCs=tf(Cs);
Ks=ss(Cs.A,Cs.A*Cs.B,Cs.C,Cs.C*Cs.B);
tfKs=tf(Ks);
%tfctrlbar=minreal(inv(1+tfctrl*Dyu)*tfctrl)%C(s)
% [Ac,Bc,Cc,Dc]=ssdata(ss(tfctrlbar));
%now compute a new controller K(s)=s*C(s)
% tf0=tf([1,0],[0,1]);
% tfk=minreal(tf0*tfctrlbar);
% sysctrlk=ss(tfk);
% syscloopk=feedback(sysoloop,sysctrlk,1)
syscloopk=feedback(sysoloop,Ks,1);
axis(gca,'normal')
pzplot(syscloopk);
[Wn,zeta] = damp(syscloopk);%coefficients satisfied!
impulseplot(sysoloop,'b--',syscloopk,'r',20)
%figure(1,1)corresponds to 'rudder to yaw rate'
%figure(2,2)corresponds to 'aileron to bank angle'
% Ak=Ac;
% Bk=A*Bc;
% Ck=Cc;
% Dk=Cc*Bc;
% sysctrlk=ss(Ak,Bk,Ck,Dk);
% tfctrlk=tf(sysctrlk);
% syscloopk=feedback(sysoloop,sysctrlk,1);
