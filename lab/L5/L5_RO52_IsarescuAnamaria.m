clear all
close all
clc
load('lab5_5.mat');

u_id = id.InputData;
y_id = id.OutputData;
u_val = val.InputData;
y_val = val.OutputData;

%plot(tid,u_id,tid,y_id);
% figure
% plot(tid,y_id);
% title("Iesire identificare");
% figure
% plot(tid,u_id);
% title("Intrare identificare");

u_id = detrend(u_id);
y_id = detrend(y_id);
figure
plot(tid,u_id);
title("Intrarea de medie zero");
figure
plot(tid,y_id);
title("Iesirea de medie zero");

l_uid = length(u_id); %N-numarul de esantioane
ru_id = [];
ryu_id = [];
for tau = 0:l_uid - 1
    %tau = tau + 1
    su = 0;
    syu = 0;
    for k = 1:l_uid - tau
        %k = k + 1;
        su = su + u_id(k + tau)*u_id(k);
        syu = syu + y_id(k + tau)*u_id(k);
    end
    ru = 1/l_uid*su;
    ru_id(end + 1) = ru;
    ryu = 1/l_uid*syu;
    ryu_id(end + 1) = ryu;
end
ryu_id = ryu_id';
% M = 10;
% T = 55;
M = 36;
T = l_uid;

Ryu_id = [];
Ru_id = [];
for i = 1:T
    Ryu_id(end + 1) = ryu_id(i);
end

for i = 1:T
    for j = 1:M
        Ru_id(i,j) = ru_id(abs(i - j) + 1);
    end
end

%calcul H_hat
H = Ru_id\Ryu_id';

% calcul iesire prezisa
y_hat_id = conv(H,u_id);

%calcul MSE pentru datele de indentificare
s_mse_id = 0;
for i = 1:length(y_id)
        s_mse_id = s_mse_id + (y_hat_id(i) - y_id(i)).^2;
end
mse_identificare = (1/length(y_id))*s_mse_id


figure
plot(H,'r');
title("H vs imp");
hold on
plot(imp,'b')


figure
plot(y_hat_id(1:length(y_id)),'r');
hold on
plot(y_id,'b');
title("y hat id vs y id");

%calcul y_hat petru validare
y_hat_val = conv(H,u_val);

%calcul MSE pentru datele de validare
s_mse_val = 0;
for i = 1:length(y_val)
        s_mse_val = s_mse_val + (y_hat_val(i) - y_val(i)).^2;
end
mse_validare = (1/length(y_val))*s_mse_val

figure
plot(y_hat_val(1:length(y_val)),'r');
hold on
plot(y_val,'b');
title("y hat val vs y val");
legend('y hat val','y val');
