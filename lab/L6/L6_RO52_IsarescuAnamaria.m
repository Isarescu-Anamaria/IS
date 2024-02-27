close all
clear all
clc

load('lab6_3 (1).mat');
u_id = id.InputData;
y_id = id.OutputData;
u_val = val.InputData;
y_val = val.OutputData;

%plot(id);
figure
subplot(211);
plot(u_id);
title('Date de intrare de identificare');
subplot(212);
plot(y_id);
title('Date de iesire de identificare');

figure
subplot(211);
plot(u_val);
title('Date de intrare de validare');
subplot(212);
plot(y_val);
title('Date de iesire de validare');

na = 10;
nb = 10;
%predictie
%identificare
phi_id = [];
N_id = length(y_id);
for i = 1:N_id
    %for pt y
    for j = 1:na 
        if i-j <= 0
            phi_id(i,j) = 0;
        else
            phi_id(i,j) = -y_id(i-j);
        end
        
    end

    %for pt u
    for j = na + 1:na + nb
        c = j - na;
        if i-c <= 0
            phi_id(i,j) = 0;
        else
            phi_id(i,j) = u_id(i-c);
        end
    end
end

theta = phi_id\y_id;
y_id_hat = phi_id*theta;

% predictie pentru identificare
figure
plot(y_id_hat);
hold on
plot(y_id);
title('Predictie pentru idetificare');

%mse 
mse_id_p = 1/N_id*sum((y_id_hat - y_id).^2);

%validare
phi_val = [];
N_val = length(y_val);
for i = 1:N_val
    %for pt y
    for j = 1:na 
        if i-j <= 0
            phi_val(i,j) = 0;
        else
            phi_val(i,j) = -y_val(i-j);
        end  
    end

    %for pt u
    for j = na + 1:na + nb
        c = j - na;
        if i-c <= 0
            phi_val(i,j) = 0;
        else
            phi_val(i,j) = u_val(i-c);
        end
    end
end

y_val_hat = phi_val*theta;

% predictie pentru validare
figure
plot(y_val_hat);
hold on
plot(y_val);
title('Predictie pentru validare');

%mse 
mse_val_p = 1/N_val*sum((y_val_hat - y_val).^2)

%simulare pentru datele de validare
y_hat_simulare = zeros(length(y_val),1);
%y_hat_simulare = [];
 for i = 1:N_val
    for j = 1:na
        if i-j > 0
            y_hat_simulare(i) = y_hat_simulare(i) - theta(j)*y_hat_simulare(i-j);
        end
    end

    for j = 1:nb
        if i-j > 0
            y_hat_simulare(i) = y_hat_simulare(i) + theta(j+na)*u_val(i-j);
        end
    end
 end

figure
plot(y_val);
hold on
plot(y_hat_simulare)
title('Simulare pentru validare');

%mse
N_val_simulare = length(y_hat_simulare);
mse_val_s = 1/N_val_simulare*sum((y_hat_simulare - y_val).^2)

