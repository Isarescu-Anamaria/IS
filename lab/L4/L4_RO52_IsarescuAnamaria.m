load('lab4_order1_3.mat');
u_ord1 = data.InputData;
y_ord1 = data.OutputData;
figure
plot(t,u_ord1,t,y_ord1);
hold on
title('Date de identificare pentru sistemul de ordinul 1');
yss_ord1 = mean(y_ord1(101:130));
uss_ord1 = u_ord1(1);
% factor de proportionalitate -> k
k = yss_ord1/uss_ord1
ymax = 9.08496;
t0 = 7.44;%timpul corespunzator lui ymax
yval = 0.368 * (ymax - yss_ord1) + yss_ord1;
plot(yval,'*b');
t1 = 10.08;
T = t1 - t0
H1 = tf(k,[T 1])
A = -1/T;
B = k/T;
C = 1;
D = 0;
Hss1 = ss(A,B,C,D);
figure
lsim(Hss1,u_ord1,t,yss_ord1);
hold on
plot(t,y_ord1,'r');

% validare

figure
lsim(Hss1,u_ord1,t,yss_ord1);
hold on
plot(t(131:330),y_ord1(131:330),'r');
%mse pentru datele de validare
t_validare = 131:330;
u_validare = u_ord1(131:330);
y_aproximat1 = lsim(Hss1,u_validare,t_validare,yss_ord1);
eroare = y_ord1(131:330) - y_aproximat1;
S = sum(eroare.^2);
mse1 = S/length(t_validare)

%%
load('lab4_order2_3.mat');
u_ord2 = data.InputData;
y_ord2 = data.OutputData;
figure
plot(t,u_ord2,t,y_ord2);
title('Date de identificare pentru sistemul de ordinul 2');
figure
plot(t,y_ord2);
hold on

% calcul k
yss_ord2 = mean(y_ord2(101:130));
uss_ord2 = u_ord2(1);
k = yss_ord2/uss_ord2

%suprareglaj

%plot(yss_ord2,'*r')
%t_yss = 0.933333;
t00 = 3.86667;
t01 = 8;
t02 = 12.4;
Ts = t(2) - t(1);
k00 = 31; %t00/Ts
k01 = 64; %t01/Ts
k02 = 90; %t02/Ts
s1 = sum(y_ord2(k00:k01) - yss_ord2);
A = Ts*s1;
s2 = sum(yss_ord2 - y_ord2(k01:k02));
a = Ts*s2;
M = a/A;
tita = log(1/M)/sqrt(pi^2 + log(M)^2)
%calcul T0
t1 = 5.2;
t3 = 12.8;
T0 = t3 - t1;
%calcum wn
wn = 2*pi/(T0*sqrt(1 - tita^2))
num = k*wn^2;
den = [1 2*tita*wn wn^2];
H2 = tf(num,den)

A = [0 1;-wn^2 -2*tita*wn];
B = [0;k*wn^2];
C = [1 0];
D = 0;
Hss2 = ss(A,B,C,D);
% validare
figure
lsim(Hss2,u_ord2,t,[yss_ord2 0]);
hold on
%plot(t,y_ord2,'r');
plot(t(131:330),y_ord2(131:330),'r');

%mse pentru datele de validare

t_validare2 = 131:330;
u_validare2 = u_ord2(131:330);
y_aproximat2 = lsim(Hss2,u_validare2,t_validare2,[yss_ord2 0]);
eroare2 = y_ord2(131:330) - y_aproximat2;
S2 = sum(eroare2.^2);
mse2 = S2/length(t_validare2)