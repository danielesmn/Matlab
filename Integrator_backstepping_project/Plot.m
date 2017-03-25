clc 
clear all
close all force

%% This code implement the optimization backward in time to estimate the actual time of the state (parameters) of a dynamical system discretized by backward Euler
tic

global YY Nwindow  time Ts interN startbacktime actualk

N = 2000;
Nwindow = 10; 
Ts = 0.005; %sampling time
time = 0:Ts:(N-1)*Ts;
interN = 30; %intersampling parameter

%INITIAL CONDITIONS
x10 = 8.3;
x20 = 6.2;
x30 = -14.5;

x10_hat = 8.3;
x20_hat = 6.2;
x30_hat = -14.5;

reconstruction_error = 1E-5;

myfontsize = 22;

xs = zeros(3,N);
xhat = zeros(3,N);
xold = zeros(3,N);
y = zeros(1,N);
yhat  = zeros(1,N);  %estimated output
e = zeros(1,N); % reconstruction error
G = zeros(1,Nwindow); %window data
J = zeros(1,ceil(N/interN)); % cost function
E = zeros(1,ceil(N/interN)); % ERROR TO EVALUATE COST FUNCTION
Y = zeros(1,ceil(N/interN)); % INTERSAMPLED OUTPUT
A = zeros(1,N);
B = zeros(1,N);
eta_z = zeros(1,N);
u = zeros(1,N);

xs(:,1) = [x10;x20;x30];
xhat(:,1) = [x10_hat;x20_hat;x30_hat];
%%

j=1; %counter of intersamples
Y(j) = y(1);
E(j) = e(1);
J(j) = E(max(1-Nwindow+1,1):1)*E(max(1-Nwindow+1,1):1)';  
j=2;
startbacktime = 1;

for k=2:N
    
    gg = 10;
    
    eta_z(k) = (2*xs(1,k-1) + 2*xs(2,k-1) + xs(3,k-1));

    A(k) = A(k-1) + Ts*(-xs(3,k-1)*eta_z(k-1));
    B(k) = B(k-1) + Ts*(-xs(2,k-1)*eta_z(k-1));

    u(k) = (-4*xs(1,k-1) - (5-B(k-1))*xs(2,k-1) - (3-A(k-1))*xs(3,k-1) + xs(1,k-1)^3 -gg*eta_z(k-1));
    
    %dinamica dell'impianto discretizzata
    xs(1,k) = xs(1,k-1) + Ts*(xs(2,k-1));
    xs(2,k) = xs(2,k-1) + Ts*(xs(3,k-1));
    xs(3,k) = xs(3,k-1) + Ts*(-0.7*xs(3,k-1) - 1*xs(2,k-1) + xs(1,k-1) - xs(1,k-1)^3 + u(k-1));
    
    %uscita dell'impianto
    y(k) =   xs(1,k) ;%change
    
    xhat(1,k) = xhat(1,k-1) + Ts*(xhat(2,k-1));
    xhat(2,k) = xhat(2,k-1) + Ts*(xhat(3,k-1));
    xhat(3,k) = xhat(3,k-1) + Ts*(-0.7*xhat(3,k-1) - 1*xhat(2,k-1) + xhat(1,k-1) - xhat(1,k-1)^3 + u(k-1));
    
    yhat(k) = xhat(1,k);
    
    e(k) =  y(k) - yhat(k);
    
    %intersample a value
    if(mod(k,interN)==1)
       
        Y(j) = y(k);
        E(j) = e(k);
        J(j) = E(max(j-Nwindow+1,1):j)*E(max(j-Nwindow+1,1):j)';  
      
        if(j >= Nwindow)
            
            if(abs(E(j))>reconstruction_error)
                
                YY =  Y(j-Nwindow+1:j);
                actualk = k;
                xhat(:,k) = fminsearch('funzione_costoJ',xhat(:,k));
                
                %Re-evaluate the estimated state (and all other
                %variables) from startbacktime to k with the new value
                %provided by fmincon
                jj = j-Nwindow;        
            end         
        end
        j = j + 1;
    end  
end         
%%
figure(1)
ax(1) = subplot(4,1,1);
plot(time,y,time,yhat,'LineWidth',1);
grid on
legend('y','yhat');
ax(2) = subplot(4,1,2);
plot(time,e,'LineWidth',1)
grid on
ylabel('$e$','Interpreter','Latex','FontSize',myfontsize);
ax(3) = subplot(4,1,3);
plot(time(1:interN:end),J,'LineWidth',1);
grid on
ylabel('$J$','Interpreter','Latex','FontSize',myfontsize);
ax(4) = subplot(4,1,4);
plot(time,xs,'LineWidth',1)
grid on
hold on
plot(time,xhat,'--','LineWidth',1)
ylabel('x','Interpreter','Latex','FontSize',myfontsize);
linkaxes(ax,'x')

%%
figure(2)
plot(time,xhat-xs)
grid on
ylabel('error','FontSize',myfontsize);

%%
figure(3)
plot(time,xhat(1, :),time,xhat(2, :),time,xhat(3, :),'LineWidth',1)
grid on
legend('x1hat','x2hat','x3hat');
ylabel('estimate x','FontSize',myfontsize);

toc

figure(4)

x0 = [0 0 0.1];
%x0_hat = [8.3, 6.2, -14.5];
%x0_hat = [1, 1, -1];
x0_hat = [8.3, 6.2, -14.5];

K = 3;

sim('Backstepping');
        
%phase_plane(x,'r');
phase_plane(ScopeData1,'b');
phase_plane1(xhat,'g');
legend('X High gain Observer','X Newton like Observer');


