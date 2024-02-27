clear all
close all
clc

u_spab = Spab(10,300,-0.8,0.8);

figure
plot(u_spab);
title("Iesire SPAB 10");

u_z = zeros(10,1);
u_step = 0.3*ones(70,1);
u = [u_z;u_spab;u_z;u_step];
figure
plot(u);
title("Date de intrare");

% [vel,alpha,t] = run(u,'10')
% load('motor9.mat');
load('date_lab9.mat');
figure
plot(vel);

na = 4;
nb = 4;
nk = 1;
% u_id = u(1:320);
% y_id = vel(1:320);
% u_val = u(321:390);
% y_val = vel(321:390);

u_id = u(1:360);
y_id = vel(1:360);
u_val = u(450:650);
y_val = vel(450:650);
Ts = 0.01;

id = iddata(y_id',u_id,Ts);
val = iddata(y_val',u_val,Ts);
model_arx = arx(id, [na, nb, nk]);
figure
compare(model_arx, id);
title("Arx pentru identificare");
y = compare(model_arx, id);
y_hat = y.OutputData;
% A = model_arx.A;
% B = model_arx.B;
N = length(y_id);
%variabile instrumentale
Z = [];
for k = 1:length(y_id)
    z = [];
    for i = 1:na
        if k - i <= 0
            z(end + 1) = 0;
        else
            z(end + 1) = - y_hat(k - i);
        end
    end
    for i = nk:nk + nb - 1
        if k - i <= 0
            z(end + 1) = 0;
        else
            z(end + 1) = u_id(k - i);
        end
    end
    Z = [Z;z];
end

PHI = [];
for k = 1:length(y_id)
    phi = [];
    for i = 1:na
        if k - i <= 0
            phi(end + 1) = 0;
        else
            phi(end + 1) = - y_id(k - i);
        end
    end
    for i = nk:nk + nb - 1
        if k - i <= 0
            phi(end + 1) = 0;
        else
            phi(end + 1) = u_id(k - i);
        end
    end
    PHI = [PHI;phi];
end

PHI = PHI';
Z = Z';

PHI_tilda = [];
s_phi = zeros(na + nb,na + nb);

for k = 1:N
    aux = zeros(na + nb,na + nb);
    aux_z = zeros(na + nb,1);
    aux_phi = zeros(na + nb,1);
    for i = 1:na + nb
        aux_z(i) = Z(i,k);
        aux_phi(i) = PHI(i,k);
    end
    aux = aux_z .* aux_phi';
    s_phi = s_phi + aux; 
end
PHI_tilda = s_phi/N;

Y_tilda = zeros(na + nb,1);

for i = 1:na + nb
    s_y = 0;
    for k = 1:N
       s_y = s_y + Z(i,k) * y_id(k); 
    end
    Y_tilda(i) = s_y/N;
end

theta = PHI_tilda\Y_tilda;

A = [];
B = [];

A(1) = 1;
for i = 2:na + 1
    A(i) = theta(i - 1); 
end

B(1) = 0;
for i = 2:nb + 1
    B(i) = theta(na + i - 1);
end

model_vi = idpoly(A,B,1,1,1,0,Ts);
figure
compare(model_vi,val);


%%
function u = Spab(m,N,a,b)

for m = 3:10
if (m == 3)
    c = [1 0 1];
end

if (m == 4)
    c = [1 0 0 1];
end

if (m == 5)
    c = [0 1 0 0 1];
end

if (m == 6)
    c = [1 0 0 0 0 1];
end

if (m == 7)
    c = [1 0 0 0 0 0 1];
end

if (m == 8)
    c = [1 1 0 0 0 0 1 1];
end

if (m == 9)
    c = [0 0 0 1 0 0 0 0 1];
end

if (m == 10) 
    c = [0 0 1 0 0 0 0 0 0 1];
end

A = zeros(m,m);
for i = 1:m
    for j = 1:m
        if i == 1
            A(i,j) = c(j);
        elseif i-j == 1
            A(i,j) = 1;
        else
            A(i,j) = 0;
        end
    end
end
end

x = zeros(m,1);
x(1) = 1;
x(4) = 1;
u = zeros(N,1);
u(1) = x(m);
for k = 2:N
    x = mod(A*x,2);
    u(k) = x(m);
end

u = a + (b - a)*u;
end