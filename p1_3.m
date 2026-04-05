% Limpiar y cerrar ventanas emergentes
clc, close all, clear;

%% Lectura del archivo csv con los datos del motor
data = readtable("data_motor.csv");
t = data.time_t_;            % Time
u = data.ex_signal_u_;       % U signal 
y = data.system_response_y_; % Y signal response

%% Metodo de Ziegler and Nichols
y_base = min(y);                 
y_estable = mean(y(end-10:end)); 

% Obtener la ecuacion de la recta tangente
t1 = 0.757575; y1 = 0.405307;
t2 = 0.808080; y2 = 0.483201;
m = (y2 - y1) / (t2 - t1);       
b = y1 - m * t1;
t_t = linspace(0.34273, 1.3153, 100);
ylin = m*t_t + b;              

% Detectar el instante real donde el escalón se activa (no duplicar retrasos)
t_escalon = t(find(u > 0, 1, 'first')); 

% Parámetros del modelo
K = (y_estable - y_base) / (max(u) - min(u));                   
t_base = (y_base - b) / m;        % Corte de tangente con la línea base
t_100 = (y_estable - b) / m;      % Corte de tangente con la línea 100%

theta = t_base - t_escalon;       % Retardo real relativo al escalón
tau = t_100 - t_base;             % Constante de tiempo ZN

% Funcion de transferencia para Ziegler and Nichols y simulación
G1 = tf(K, [tau 1], 'InputDelay', max(0, theta));
y_zn = lsim(G1, u, t);            % Usamos 'u' directamente

%% Metodo Analitico
y_632 = y_base + 0.632 * (y_estable - y_base);
y_284 = y_base + 0.284 * (y_estable - y_base);

[~, idx_632] = min(abs(y - y_632));
[~, idx_284] = min(abs(y - y_284));
t_632 = t(idx_632); 
t_284 = t(idx_284);

tau2 = 3 * (t_632 - t_284) / 2;
theta2 = (t_632 - tau2) - t_escalon; % Ajuste relativo al escalón

G2 = tf(K, [tau2 1], 'InputDelay', max(0, theta2));
y_analitico = lsim(G2, u, t);

%% Metodo de Miller
tau3 = t_632 - t_base; % Miller usa t_632 menos el t_base de la tangente
G3 = tf(K, [tau3 1], 'InputDelay', max(0, theta));
y_m = lsim(G3, u, t);

%% FIT
fit_zn = 100 * (1 - norm(y - y_zn) / norm(y - mean(y)));
fit_analitico = 100 * (1 - norm(y - y_analitico) / norm(y - mean(y)));
fit_m = 100 * (1 - norm(y - y_m) / norm(y - mean(y)));

%% Mostrar en pantalla
fprintf('\n===============================================\n');
fprintf('        PARÁMETROS IDENTIFICADOS DEL SISTEMA\n');
fprintf('===============================================\n');
fprintf('Ganancia del proceso (K): %.4f\n\n', K);
fprintf('Método\t\t\t\tθ [s]\t\tτ [s]\n');
fprintf('-----------------------------------------------\n');
fprintf('Ziegler & Nichols\t%.5f\t\t%.5f\n', theta, tau);
fprintf('Analítico\t\t\t%.5f\t\t%.5f\n', theta2, tau2);
fprintf('Miller\t\t\t\t%.5f\t\t%.5f\n', theta, tau3);
fprintf('-----------------------------------------------\n\n');

fprintf('\n===============================================\n');
fprintf('       PORCENTAJE DE AJUSTE (FIT)\n');
fprintf('===============================================\n');
fprintf('FIT Ziegler & Nichols: %.2f%%\n', fit_zn);
fprintf('FIT Analítico: %.2f%%\n', fit_analitico);
fprintf('FIT Miller: %.2f%%\n', fit_m);

%% Grafica
figure; hold on;
plot(t, u, '-g', t, y, '-b', t, y_zn, '-y', t, y_analitico, '-r', 'LineWidth', 2.5);   
plot(t, y_m, 'Color', [0.5 0 0.5], 'LineWidth', 2.5);  

plot(t, y_base * ones(size(t)), '--r', t, y_estable * ones(size(t)), '--k', t_t, ylin, '-.m', 'LineWidth', 2);      
plot(t, y_632 * ones(size(t)), 'c--', 'LineWidth', 1);
plot(t, y_284 * ones(size(t)), '--', 'Color', [0.6 0.3 0], 'LineWidth', 1);

xlabel('Tiempo (s)'); ylabel('Amplitud');
legend('Señal Escalon', 'Respuesta del Proceso', 'Ziegler and Nichols', 'Analitico', 'Miller', 'Base: 0%', 'Sistema Estable: 100%', 'Recta Tg', '63.2%', '28.4%', 'Location', 'best');
xlim([0 5]); ylim([-0.1 1.6]); grid on;