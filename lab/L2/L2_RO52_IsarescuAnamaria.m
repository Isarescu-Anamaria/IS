load('lab2_12.mat');

%pasul 1 
figure
plot(id.X,id.Y); title('Date de identificare si date de validare');
hold on
plot(val.X,val.Y);

%pasul 3
MSE_ID = []; %vector pentru mse de identificare
MSE_VAL = []; %vector pentru mse de validare
for var = 1:25
n = var;
fi_id = [];
for i = 1:length(id.X)
    for j = 1:n
        fi_id(i,j) = id.X(i)^(j-1); %calcul matrice de regresori pentru datele de identificare
    end
end
theta_id = fi_id\id.Y'; % gasire matrice theta(cu datele de la id)

%pasul 4
y_id = fi_id*theta_id; %y cu caciula pentru id
%pasul 5
sid = 0; % suma pentru mse pentru datele de identificare
for i = 1:length(id.Y)
    sid = sid + (id.Y(i) - y_id(i)).^2;
end
N1 = length(id.X); % numarul de esantioane pentru datele de identificare
mse_id = (1/N1)*sid; % calcul formula mse pentru datele de id 
MSE_ID(var) = mse_id; % formare vector pentru mse de id

fi_val = [];
for i = 1:length(val.X)
    for j = 1:n
        fi_val(i,j) = val.X(i)^(j-1); % calcul matrice de regresori pentru datele de validare
    end
end
%theta_val = fi_val\val.Y';
%y_val = fi_val*theta_val; % y cu caciula pentru val
y_val = fi_val*theta_id; % y cu caciula pentru val
sval = 0; % suma mse pentru datele de validare
for i = 1:length(val.Y)
    sval = sval + (val.Y(i) - y_val(i)).^2;
end
N2 = length(val.X); % numar de esantioane pentru datele de validare
mse_val = (1/N2)*sval; % calcul formula mse pentru datele de validare
MSE_VAL(var) = mse_val; % formare vector mse pentru datele de validare
end

[mse_val_minim,index_val_minim] = min(MSE_VAL); %gasire valoare mse minima in vectorul de mse-uri si pozitia lui in vector
figure
plot(MSE_VAL); title('MSE pentru datele de validare'); % grafic pentru vectorul de mse-uri pentru validare
hold on
%mse_val_minim = min(MSE_VAL);
plot(index_val_minim,mse_val_minim,'*g');% afisare punct de minim
%[mse_id_minim,index_id_minim] = min(MSE_ID);
%mse_id_minim = min(MSE_ID);
%plot(mse_id_minim,'*y');

%pasul 6
figure
plot(val.X,val.Y,val.X,y_val); title('Functia aproximata pentru datele de validare');% grafic pentru y de validare vs y cu caciula de validare
