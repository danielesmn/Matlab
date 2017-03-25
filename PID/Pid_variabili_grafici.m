clc

clear all
close all

epsilon = 0.00001;
Ta_5 = 0.0504;
Tau = Ta_5/3;
T = 0.04;
Td = T/2;

zita1 = 0.46;
omega_n1 = 2 * pi/0.042;
ALF1 = tf([omega_n1^2],[1, 2 * zita1 * omega_n1, omega_n1^2]);

zita2 = 0.8;
omega_n2 = 2 * pi/0.05;
my_DEP2 = 0.7 * tf([omega_n2^2],[1, 2 * zita2 * omega_n2, omega_n2^2]);%% perchè nonostante 0.7 si assesta sempre a 1

mytime = [0:0.001:1];
input = sin(15 * mytime);

L = ALF1 * my_DEP2;

Kp_signed = 2.26;
Kp = 0.6*Kp_signed;
Kd = Kp*Td;
Ki = Kp*0.1*140;

figure(1)
lsim(L,input,mytime)
rlocus(L)
grid on

figure(2)
bode(L)
grid on

figure(3)
step(L)
grid on

figure(4)
nyquist(L)
grid on

break

