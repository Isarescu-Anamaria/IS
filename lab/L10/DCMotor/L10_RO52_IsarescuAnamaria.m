clear all
close all
clc

u_z = zeros(20,1);
u_step = 0.3*ones(70,1);

u_val = [u_z;u_step;u_z];
u_id = idinput(200,'prbs',[],[-0.8 0.8]);
na = 6;
nb = 6;
Ts = 0.01;

%%
motor = DCMRun.start("Ts", 10e-3);

for k = 1:length(u_val)
    motor.wait;
    y_val(k)= motor.step(u_val(k));
end

P_inv = [];
P_inv = 1000*eye(na + nb,na + nb);
theta_hat = zeros(na + nb,1);
eroare = [];
w = zeros(na + nb);
y_id = [];

for k = 1:length(u_id)
    motor.wait;
    y_id(k) = motor.step(u_id(k));
   
    phi = zeros(na+nb,1);
    for j = 1:na 
        if k-j <= 0
            phi(j) = 0;
        else
            phi(j) = -y_id(k-j);
        end
        
    end

    for j = na + 1:na + nb
        c = j - na;
        if k-c <= 0
            phi(j) = 0;
        else
            phi(j) = u_id(k-c);
        end
    end

    if k == 1
        eroare = y_id(k);
    else
    eroare = y_id(k) - phi'*theta_hat;
    end
    
    phi = phi';
    P_inv = P_inv - (P_inv*(phi')*phi*P_inv)/(1+phi*P_inv*(phi'));
    w = P_inv*phi';  %vector coloana de na+nb
    theta_hat = theta_hat + w*eroare;

    motor.wait();
end
motor.stop();

A = [];
B = [];

A(1) = 1;
for i = 2:na + 1
    A(i) = theta_hat(i - 1); 
end

B(1) = 0;
for i = 2:nb + 1
    B(i) = theta_hat(na + i - 1);
end

model_complet = idpoly(A,B,[],[],[],0,Ts);
val = iddata(y_val',u_val,Ts);


figure
compare(model_complet,val)

