clc
clear all
syms Rprimepara A0para Spara rp1 rp2 s0 s1 s2 a01 a02 x T r0
%% Tf Generation

alpha=0.003725;
beta=0.0005685;
c=1.5;
L0=1;
m=0.2;
k=0.2;
Lp=600*pi*(26+38.5);
R=16.2*Lp/1000;
x0=-0.001;

i0 = sqrt((-k/beta)*((x0^3)+2*alpha*(x0^2)+(alpha^2)*x0));
e = R*i0;

A=[0 1 0;
   -(k/m)-beta*(i0^2)/(m*2*((x0+alpha)^2)*x0)  -c/m  -beta*i0/(m*2*((x0+alpha)^2));
    e/beta-R*i0/(2*beta) i0/(3*((x0+alpha)))  -(R*x0/(2*beta))-(alpha*R/beta)];
B=[0;0;alpha/beta];
C=[1 0 0];
D=[0];
[b,a]=ss2tf(A,B,C,D);
sys=tf(b,a)
poles=roots(a)
pole1=double(poles(2));
pole2=double(poles(3));
syss=tf(b,[1 -(pole1+pole2) pole1*pole2])

%% Controller design

zeta=0.7;
TS=0.1;
wn=4/(zeta*TS);
sysde=tf([wn^2],[1 2*zeta*wn wn^2])


[B,A]=tfdata(syss);
A2=cell2mat(A);
B2=cell2mat(B);

% TF desired

[Bm,Am]=tfdata(sysde,'v');

Am2=poly2sym(Am);
Bm2=poly2sym(Bm);

% STR degree

degreeR=1;
degreeT=1;
degreeS=1;

degreeBminus=0;
degreeBplus=0;
degreeBmprime=0;
degreeRprime=1;
degreeA0=1;

%% STR parameters

% Bplus

Bplus=1;
Bplus2=poly2sym([1]);

% Bminus

Bminus=tf(B2(3),1);
[forsizeBminusnum,forsizeBminusden]=tfdata(Bminus);
Bminus2=B2(3);

% Bm

BM=tf(Bm(3),1);
forsizeBmnum=Bm(3);

% Bmprime

Bmprime=BM/Bminus2;
Bmprimevec=forsizeBmnum/Bminus2;

% Rprime

rprime=[1 r0];
%Rprime=tf(rprime,1);
Rprime2=poly2sym(rprime);

% A0

A0=[1 1000];
A02=poly2sym(A0);

% S

S=[s1 s0];
S2=poly2sym(S);

% diophantine

S3=Bminus2*S2;
S4=vpa(S3);

R3=poly2sym(A2)*Rprime2;
R4=vpa(R3);

A3=A02*Am2;
A4=vpa(A3);

dio= R4+S4-A4;
diop=vpa(dio);
dioph=coeffs(diop,x);


eqns=[dioph(1)==0,dioph(2)==0,dioph(3)==0];
diophantine1 = solve(eqns, [s1 s0 r0]);

Sdio0 = double(diophantine1.s0)
Sdio1 = double(diophantine1.s1)
rdio0 = double(diophantine1.r0)

A0dio=[1];
Sdio=[Sdio1 Sdio0]
Rprimedio=[1 rdio0];

ttp=(Bmprimevec*A0);
Tp=tf(ttp,1)
Sp=tf(Sdio,1)
Rpo=tf(Bplus*Rprimedio,1)

Closedloop=syss*Tp/(Rpo+Sp*syss);
fb=bandwidth(Closedloop);
Ts=0.005*2*pi/fb;

sysded=c2d(sysde,Ts);
sysd=c2d(syss,Ts);

SR=Sp/Rpo
SRd=c2d(SR,Ts,'zoh')
[sdd,rdd]=tfdata(SRd);
TR=Tp/Rpo
TRd=c2d(TR,Ts,'zoh')
[tdd,rdd]=tfdata(TRd);
sd=cell2mat(sdd)
rd=cell2mat(rdd)
td=cell2mat(tdd)

