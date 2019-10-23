clear; close all; clc
N0 = -206;
Fd = 5.2e6; % Частота дискретизации, МГц
T = 3.15e-3; % Интервал накопления, мс
Pf = 1e-5; % Вероятность ЛТ
M = 1; % Ячеек по частоте
N = 52; % Ячеек по задержке
delta_tau = 1/2; % Шаг ячеек по задержке, символов
tau_tilda = (0:N-1)*delta_tau + delta_tau/2; % Опорные задержки
tau_max = N*delta_tau;
% Число экспериментов для поиска порога 
J1 = 100000;
% Число экспериментов для расчета характеристик обнаружения
J2 = 10000;
% Определение inline-функции ro(dtau)
ro = inline('(1 - abs(dtau)) .* (abs(dtau)<1)', 'dtau');
qcno_dB = 25:0.5:45;
figure('Name','Xарактеристики обнаружения');
Td = 1/Fd;
std_y = 10^(N0/10)/(2*Td); % СКО шума выборки
L = T / Td; % Число суммирований
std_IQ = std_y * sqrt(L/2); % СКО шума корреляционных сумм
delta_f = 2/3 / T; % Шаг ячеек по частоте, Гц
omega_tilda = 2*pi*((0:M-1)*delta_f + delta_f/2); % Опорные частоты
omega_max = 2*pi*M*delta_f;
%Сетка, задающаяцентрыячеек
[tau_tilda_m, omega_tilda_m] = meshgrid(tau_tilda, omega_tilda);
X2max = nan(1, J1); % Инициализация памяти
for j = 1:J1
    I = std_IQ*randn(M, N) ;
    Q = std_IQ*randn(M, N) ;
    X2 = I.^2 + Q.^2;
    X2max(j) = max(max(X2));
end
R = std_IQ^2; % Очень низкий порог
while sum(X2max > R) / J1 > Pf
    R = R * 1.0005; % Увеличиваем на 0.002 дБ
end
X2max = nan(1, J2); % Стирание прошлых результатов
Pd = nan(1, length(qcno_dB));
for q = 1:length(qcno_dB)
    qcno = 10^(qcno_dB(q)/10); % Перевод из дБ в разы
    A = 2*std_y * sqrt(qcno*Td); % Расчет амплитуды для данного с/ш
    % Истинная задержка для каждого эксперимента
    tau = tau_max * rand(1, J2);
    % Истинная частота для каждого эксперимента
    omega = omega_max * rand(1, J2);
    % Начальная фаза в каждом эксперименте случайна
    dphi = 2*pi*rand(1, J2);
    for j = 1:J2
        dtau = tau(j) - tau_tilda_m;
        domega = omega(j) - omega_tilda_m;
        I = A*L/2 * ro(dtau) .* sinc(domega*T/2 /pi) .* cos(domega*T/2+ dphi(j)) + std_IQ*randn(M, N) ;
        Q =- A*L/2 * ro(dtau) .* sinc(domega*T/2 /pi) .* sin(domega*T/2+ dphi(j)) + std_IQ*randn(M, N) ;
        X2 = I.^2 + Q.^2;
        X2max(j) = max(max(X2));
    end
    Pd(q) = sum(X2max > R) / J2; % Среднее превышение порога в
end
plot(qcno_dB, Pd);
xlabel('q_{c/n0}, dBHz')
ylabel('P_d        ','Rotation',0)
grid on