clear all
close all
clc

u_spab = Spab(10,200,-0.7,0.7);

figure
plot(u_spab);
title("Iesire SPAB 10");

u_z = zeros(10,1);
%u_step = 0.4*ustep(70);
u_step = 0.4*ones(70,1);
u = [u_z;u_spab;u_z;u_step];
figure
plot(u);
title("Date de intrare");

%[vel,alpha,t] = run(u,'6');
load('motor8.mat');
figure
plot(vel)
title("Date de iesire");

y = vel;
N = length(vel);
e = zeros(1,N);
d_theta = zeros(2,N);
nk = 2;

f = 1;
b = 1;
pas = 0.1; %alfa
prag = 0.001;%pragul de convergenta(delta)
cmax = 38; 
c = 0; %counter
theta_anterior = [1;1];
theta = [1;1];
vector_norme = [];

for k = 1:nk
    e(k) = y(k);
    d_theta(1,k) = 0;
    d_theta(2,k) = 0;
end

for k = nk + 1:N
    e(k) = y(k) + f * y(k - 1) - b * U(k - nk) - f * e(k - 1);
    d_theta(1,k) = y(k - 1) - e(k - 1) - f * d_theta(1,k - 1);
    d_theta(2,k) = -U(k - nk) - f * d_theta(2,k - 1);
end

%gradientul functiei V
s1_g = 0;
s2_g = 0;
for k = 1:N
    s1_g = s1_g + e(k)*d_theta(1,k);
    s2_g = s2_g + e(k)*d_theta(2,k);
end
dV_theta = zeros(2,1);
dV_theta(1,1) = (2/(N - nk))*s1_g;
dV_theta(2,1) = (2/(N - nk))*s2_g;

%%Hessianul
H = zeros(2,2);
h = zeros(2,1);
for k = 1:N
    h = [d_theta(1,k);d_theta(2,k)];
    H = H + h * h';
end
H = (2/(N - nk))*H;

%formula de actualizare Gauss-Newton
theta = theta_anterior - pas * inv(H) * dV_theta;
c = c + 1;
norma = norm(theta - theta_anterior);
vector_norme(end+1) = norma;


while norma <= prag || c < cmax
e = zeros(1,N);
d_theta = zeros(2,N);
f = theta(1,1);
b = theta(2,1);
theta_anterior = theta;
for k = 1:nk
    e(k) = y(k);
    d_theta(1,k) = 0;
    d_theta(2,k) = 0;
end

for k = nk + 1:N
    e(k) = y(k) + f * y(k - 1) - b * U(k - nk) - f * e(k - 1);
    d_theta(1,k) = y(k - 1) - e(k - 1) - f * d_theta(1,k - 1);
    d_theta(2,k) = -U(k - nk) - f * d_theta(2,k - 1);
end

%gradientul functiei V
s1_g = 0;
s2_g = 0;
for k = 1:N
    s1_g = s1_g + e(k)*d_theta(1,k);
    s2_g = s2_g + e(k)*d_theta(2,k);
end
dV_theta = zeros(2,1);
dV_theta(1,1) = 2/(N - nk)*s1_g;
dV_theta(2,1) = 2/(N - nk)*s2_g;

%%Hessianul
H = zeros(2,2);
h = zeros(2,1);

for k = 1:N
    h = [d_theta(1,k);d_theta(2,k)];
    H = H + h * h';
end

H = (2/(N - nk))*H;

%formula de actualizare Gauss-Newton
theta = theta_anterior - pas * inv(H) * dV_theta;
c = c + 1;
norma = norm(theta - theta_anterior);
vector_norme(end+1) = norma;
end
f = theta(1,1);
b = theta(2,1);
Ts = 0.01;
val = iddata(vel(231:300)',U(231:300)',Ts);
model = idpoly(1,[zeros(1,nk),b],1,1,[1 f],0,Ts); 
figure
compare(val,model);
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
x(4) = 1;
u = zeros(N,1);
u(1) = x(m);
for k = 2:N
    x = mod(A*x,2);
    u(k) = x(m);
end

u = a + (b - a)*u;
end