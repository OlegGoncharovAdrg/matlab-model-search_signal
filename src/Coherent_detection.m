clear; close all; clc
N0 = -206;
Fd = 5.2e6; % ������� �������������, ���
T = 3.15e-3; % �������� ����������, ��
Pf = 1e-5; % ����������� ��
M = 1; % ����� �� �������
N = 52; % ����� �� ��������
delta_tau = 1/2; % ��� ����� �� ��������, ��������
tau_tilda = (0:N-1)*delta_tau + delta_tau/2; % ������� ��������
tau_max = N*delta_tau;
% ����� ������������� ��� ������ ������ 
J1 = 100000;
% ����� ������������� ��� ������� ������������� �����������
J2 = 10000;
% ����������� inline-������� ro(dtau)
ro = inline('(1 - abs(dtau)) .* (abs(dtau)<1)', 'dtau');
qcno_dB = 25:0.5:45;
figure('Name','X������������� �����������');
Td = 1/Fd;
std_y = 10^(N0/10)/(2*Td); % ��� ���� �������
L = T / Td; % ����� ������������
std_IQ = std_y * sqrt(L/2); % ��� ���� �������������� ����
delta_f = 2/3 / T; % ��� ����� �� �������, ��
omega_tilda = 2*pi*((0:M-1)*delta_f + delta_f/2); % ������� �������
omega_max = 2*pi*M*delta_f;
%�����, �������������������
[tau_tilda_m, omega_tilda_m] = meshgrid(tau_tilda, omega_tilda);
X2max = nan(1, J1); % ������������� ������
for j = 1:J1
    I = std_IQ*randn(M, N) ;
    Q = std_IQ*randn(M, N) ;
    X2 = I.^2 + Q.^2;
    X2max(j) = max(max(X2));
end
R = std_IQ^2; % ����� ������ �����
while sum(X2max > R) / J1 > Pf
    R = R * 1.0005; % ����������� �� 0.002 ��
end
X2max = nan(1, J2); % �������� ������� �����������
Pd = nan(1, length(qcno_dB));
for q = 1:length(qcno_dB)
    qcno = 10^(qcno_dB(q)/10); % ������� �� �� � ����
    A = 2*std_y * sqrt(qcno*Td); % ������ ��������� ��� ������� �/�
    % �������� �������� ��� ������� ������������
    tau = tau_max * rand(1, J2);
    % �������� ������� ��� ������� ������������
    omega = omega_max * rand(1, J2);
    % ��������� ���� � ������ ������������ ��������
    dphi = 2*pi*rand(1, J2);
    for j = 1:J2
        dtau = tau(j) - tau_tilda_m;
        domega = omega(j) - omega_tilda_m;
        I = A*L/2 * ro(dtau) .* sinc(domega*T/2 /pi) .* cos(domega*T/2+ dphi(j)) + std_IQ*randn(M, N) ;
        Q =- A*L/2 * ro(dtau) .* sinc(domega*T/2 /pi) .* sin(domega*T/2+ dphi(j)) + std_IQ*randn(M, N) ;
        X2 = I.^2 + Q.^2;
        X2max(j) = max(max(X2));
    end
    Pd(q) = sum(X2max > R) / J2; % ������� ���������� ������ �
end
plot(qcno_dB, Pd);
xlabel('q_{c/n0}, dBHz')
ylabel('P_d        ','Rotation',0)
grid on