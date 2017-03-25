clc
clear all
close all

% VALORI COSTANTI ELETTRICHE E MECCANICHE DEL NATANTE

Ra = 0.30;         %Ohm
La = 0.025;        %Henry
Km = 1;            %N*m/A
m = 1*10^5;        %Kg
g = 9.81;          %m/s^2
J = 0.16*10^7;     %Kg*m^2
f = 0.9*10^5;      %N*m*s/rad
a = 3;             %m
b = 3.74;          %m
h = m*g*(b-a);



F1 = zpk([],[-12 -0.0281+0.6727*1i -0.0281-0.6727*1i],1) % Funzione di trasferimento tra tensione di ingresso del motore(ingresso) e angolo di roll(uscita)
         
figure(1)
hold on
margin(F1);
grid on
legend('Funzione di trasferimento del sistema');
hold off

figure(2)
hold on
nyquist(F1);
grid on
legend('Funzione di trasferimento del sistema');
hold off

F2 = zpk([],[-0.0281+0.6727*1i -0.0281-0.6727*1i],1)    % Funzione di trasferimento tra coppia applicata dal moto ondoso(disturbo) e angolo di roll(uscita)

C = zpk([-0.2-0.6*i -0.2+0.6*i],[-450 -450],10^7)      % Controllore

L = C*F1;                                               % Funzione ad anello aperto
figure(3)
hold on
margin(L);
grid on
legend('Funzione ad anello aperto');
hold off

S = feedback(1,L);                                      % Funzione di sensitività
W = feedback(L,1);                                      % Funzione di trasferimento ingresso-uscita con controllore in retroazione

figure(4)
bode(L);
hold on
bode(S);
grid on
legend('Funzione a ciclo Aperto','Funzione di sensitività');
hold off
                                              
figure(5)
hold on
t1 = 0:0.1:70;
impulse(140*W,t1);                                      % Angolo di roll in uscita dal sistema con controllore in retroazione con ingresso impulsivo e disturbo nullo
grid on
legend('Risposta del sistema controllato con ingresso impulsivo e disturbo nullo');
hold off
                                          
figure(6)                                               % Angolo di roll in uscita dal sistema con controllore in retroazione con ingresso nullo e disturbo sinusoidale
hold on
t = 0:0.1:250;
Alpha = (F2/(F1*C))*tf(1,[1 15000]);
Alpha1= Alpha+W
input = 0.0001*sin(0.628*t);
[y,t] =  lsim(Alpha1*180/pi,input,t);
plot(t,y)
grid on
legend('Risposta del sistema controllato con ingresso nullo e disturbo sinusoidale');
hold off




% VARIAZIONI PARAMETRICHE KARITONOV DEL 10%

P = [1 12.0562 1.1277 5.4396];


P1 = [1*1.10 12.0562*1.10 1.1277*0.9 5.4396*0.9];
P2 = [1*0.9 12.0562*0.9 1.1277*1.10 5.4396*1.10];
P3 = [1*1.10 12.0562*0.9 1.1277*0.9 5.4396*1.10];
P4 = [1*0.9 12.0562*1.10 1.1277*1.10 5.4396*0.9];

roots(P1)
roots(P2)
roots(P3)
roots(P4)

% VARIAZIONI PARAMETRICHE KARITONOV DEL 50%

Q1 = [1*1.50 12.0562*1.50 1.1277*0.5 5.4396*0.5]
Q2 = [1*0.5 12.0562*0.5 1.1277*1.50 5.4396*1.50]
Q3 = [1*1.50 12.0562*0.5 1.1277*0.5 5.4396*1.50]
Q4 = [1*0.5 12.0562*1.50 1.1277*1.50 5.4396*0.5]

roots(Q1)
roots(Q2)
roots(Q3)
roots(Q4)

figure(9)
nyquist(L)
grid on

