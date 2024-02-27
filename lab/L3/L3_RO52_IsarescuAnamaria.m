% sistem de ordinul 1
load('lab3_order1_2.mat');
u_ord1 = data.InputData;
y_ord1 = data.OutputData;
figure
plot(t,u_ord1,t,y_ord1);
hold on
title('Sistem de ordinul 1');

%identificare

y0_ord1 = y_ord1(1);
yss_ord1 = mean(y_ord1(60:100));%media pe 41 esnatioane pentru yss 
u0_ord1 = 0;
uss_ord1 = u_ord1(1);
k = (yss_ord1 - y0_ord1)/(uss_ord1 - u0_ord1);%calcul factor de proportionalitate
val_y = y0_ord1 + 0.632*(yss_ord1 - y0_ord1); % valoarea de pe grafic cand y este 0.632 din iesire
plot(val_y,'*b');%afisare punct
t0_ord1 = 0;
t1_ord1 = 0.91;%timpul celui mai apropiat punct de pe grafic fata de val_y
T_ord1 = t1_ord1 - t0_ord1; %calcul consanta de timp a polului
%T_ord1 = 5;
H1 = tf(k,[T_ord1,1]);%gasire functie de transfer pentru sistemul de ordinul 1
figure
lsim(H1,u_ord1,t);
figure
u_id = u_ord1(1:100);%declarare date de intrare pentru prima treapta
t_id = t(1:100);%declarare timp de simulare pentru prima treapta
lsim(H1,u_id,t_id);% simulare pentru prima treapta
hold on
% afisare valori intermediare k si T
plot(k,'*g');
hold on
plot(T_ord1,'*b');

%validare

figure
plot(t,y_ord1,'r');
hold on
lsim(H1,u_ord1,t);
u_validare = u_ord1(201:500);%declarare date de validare-treptele[3:5]
t_validare = t(201:500);%declarare timp de validare
figure
plot(t_validare,y_ord1(201:500),'r');
hold on
lsim(H1,u_validare,t_validare);%simulare pentru datele de validare
%MSE pemtru datele de validare
y_aproximat = lsim(H1,u_validare,t_validare);
eroare = y_ord1(201:500) - y_aproximat;
Sigma = sum(eroare.^2);
mse1 = Sigma/length(t_validare)

%%
% sistem de ordinul 2
load('lab3_order2_2.mat');
u_ord2 = data.InputData;
y_ord2 = data.OutputData;
figure
plot(t,u_ord2,t,y_ord2);
title('Sistem de ordinul 2');

% identificare

y0_ord2 = y_ord2(1);
yss_ord2 = mean(y_ord2(80:100));%media pentru yss 
u0_ord2 = 0;
uss_ord2 = u_ord2(1);
k = (yss_ord2 - y0_ord2)/(uss_ord2 - u0_ord2);
%suprareglaj -> M
y_max = max(y_ord2(1:100)); %y(t1)
M = (y_max - yss_ord2)/(yss_ord2 - y0_ord2);
zeta = log(1/M)/sqrt(pi^2 + log(M)^2);%determinare factor de amortizare
%pulsatie naturala -> wn
t1_ord2 = 3.26667;
y_max2 = 6.28744;
t3_ord2 = 9.33333;
T0 = t3_ord2 - t1_ord2;
wn = 2*pi/(T0*sqrt(1 - zeta^2));
%gasire functie de transfer pentru sistemul de ordinul 2
num = k*wn^2;
den = [1 2*zeta*wn wn^2];
H2 = tf(num,den);

figure
lsim(H2,u_ord2,t);
figure
u_id2 = u_ord2(1:100);%declarare date de intrare pentru prima treapta
t_id2 = t(1:100);%declarare timp de simulare pentru prima treapta
lsim(H2,u_id2,t_id2);% simulare pentru prima treapta
hold on
% afisare valori intermediare M si T0
plot(M,'*g');
hold on
plot(T0,'*b');

%validare

figure
plot(t,y_ord2,'r');
hold on
lsim(H2,u_ord2,t);
u_validare2 = u_ord2(201:500);%declarare date de validare-treptele[3:5]
t_validare2 = t(201:500);%declarare timp de validare
figure
plot(t_validare2,y_ord2(201:500),'r');
hold on
lsim(H2,u_validare2,t_validare2);%simulare pentru datele de validare
%MSE pemtru datele de validare
y_aproximat2 = lsim(H2,u_validare2,t_validare2);
eroare2 = y_ord2(201:500) - y_aproximat2;
Sigma2 = sum(eroare2.^2);
mse2 = Sigma2/length(t_validare2)



