%parameters
clear;clc;
global th_max delta_e_max ath1 ath2 ath3 omega_n_th kp_th zeta_th kd_th Kdc_th ah1 Wh omega_n_h ki_h zeta_h kp_h
% Pitch Loop
th_max = 10 * pi / 180;
delta_e_max = 30 * pi / 180;

ath1 = 0.668;
ath2 = 1.27;
ath3 = -2.08;

kp_th = -delta_e_max / th_max
omega_n_th = sqrt(ath2 + kp_th * ath3)
zeta_th = 0.7
kd_th = (2 * zeta_th * omega_n_th - ath1) / ath3
Kdc_th = kp_th * ath3 / (ath2 + kp_th * ath3)

% Altitude Loop
ah1 = 830;
Wh = 20
omega_n_h = 1 / Wh * omega_n_th
ki_h = omega_n_h ^ 2 / (Kdc_th * ah1)
zeta_h = 0.7
kp_h = 2 * zeta_h * omega_n_h / (Kdc_th * ah1)
