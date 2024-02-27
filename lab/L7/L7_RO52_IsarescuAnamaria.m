clear all
close all
clc

u_spab1 = Spab(3,200,-0.7,0.7);
u_spab2 = Spab(10,200,-0.7,0.7);

figure
plot(u_spab1);
title("Iesire SPAB 3");

figure
plot(u_spab2);
title("Iesire SPAB 10");

u_z = zeros(10,1);
u_step = 0.4*ustep(70);
u = [u_z;u_spab1;u_z;u_spab2;u_z;u_step];
figure
plot(u);
title("Date de intrare");
%%
% [vel, alpha, t] = run(u, '3');

load('data.mat');

figure
plot(vel)
title("Date de iesire")

u_id3 = u(11:210);
u_id10 = u(221:420);
u_val = u(431:510);

y_id3 = vel(11:210);
y_id10 = vel(221:420);
y_val = vel(431:510);

%%
%arx pentru m = 3
% Na_val = 30;
% Nb_val = 30;
% Nk_val = 1;
Na_3 = 66;
Nb_3 = 66;
Nk_3 = 1;
Ts = 0.01;
id_3 = iddata(y_id3',u_id3,Ts);
val = iddata(y_val',u_val,Ts);
model3_identificare = arx(id_3, [Na_3, Nb_3, Nk_3]);
% model_validare = arx(val,[Na_val, Nb_val, Nk_val]);
figure
compare(model3_identificare, val);
title("Arx pentru identificare pentru m = 3");
% figure
% compare(model_validare,val);
% title("Arx pentru validare");
%%
%arx pt m = 10
Na_10 = 69;
Nb_10 = 69;
Nk_10 = 1;
id_10 = iddata(y_id10',u_id10,Ts);
model10_identificare = arx(id_10, [Na_10, Nb_10, Nk_10]);
figure
compare(model10_identificare, val);
title("Arx pentru identificare pentru m = 10");

PE_3 = 2^3 - 1;
PE_10 = 2^10 -1;
%%
function u = Spab(m,N,a,b)

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


x = zeros(m,1);
x(1) = 1;
x(2) = 1;
u = zeros(N,1);
u(1) = x(m);
for k = 2:N
    x = mod((A*x),2);
    u(k) = x(m);
end

u = a + (b - a)*u;
end
