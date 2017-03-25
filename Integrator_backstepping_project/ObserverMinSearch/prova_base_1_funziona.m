%% This code implement the optimization backward in time to estimate the actual time of the state (parameters) of a dynamical system discretized by backward Euler

clc
clear all
close all

tic

N = 1000;
Nwindow = 10; 
Ts = 0.005; %sampling time
time = 0:Ts:(N-1)*Ts;
interN = 30; %intersampling parameter

%INITIAL CONDITIONS
x10 = 0.3;
x20 = 2;

reconstruction_error = 1E-5;

lb = [-100, -100];
ub = [100, 100];

options = optimset('LargeScale','off',...
                   'MaxIter',1000,...
                   'MaxFunEvals',1000,'display','off','Algorithm','active-set');


myfontsize = 24;


x = zeros(2,N);
xhat = zeros(2,N);
xold = zeros(2,N);
y = zeros(1,N);
yhat  = zeros(1,N);  %estimated output
e = zeros(1,N); % reconstruction error
G = zeros(1,Nwindow); %window data
J = zeros(1,ceil(N/interN)); % cost function
E = zeros(1,ceil(N/interN)); % ERROR TO EVALUATE COST FUNCTION
Y = zeros(1,ceil(N/interN)); % INTERSAMPLED OUTPUT



x (:,1) = [x10;x20];
xhat (:,1) = 0*[x10;x20];

%vettore di ingressi
u = 1.5*sin(1*2*pi*time);

%%

global YY Nwindow  time Ts interN startbacktime actualk

j=1; %counter of intersamples
Y(j) = y(1);
E(j) = e(1);
J(j) = E(max(1-Nwindow+1,1):1)*E(max(1-Nwindow+1,1):1)';  
j=2;
startbacktime = 1;

for k=2:N
   
    %dinamica dell'impianto discretizzata
    x(:,k) = x(:,k-1) + Ts*( [-x(2,k-1);x(1,k-1)] ); %change   
    %uscita dell'impianto
    y(k) =   x(1,k) ;%change
    
    xhat(:,k) = xhat(:,k-1) + Ts*( [-xhat(2,k-1);xhat(1,k-1)] ); %change   
    yhat(k) = xhat(1,k); %change
   
    e(k) =  y(k)-yhat(k);
    
    %intersample a value
    if(mod(k,interN)==1)
       
        Y(j) = y(k);
        E(j) = e(k);
        J(j) = E(max(j-Nwindow+1,1):j)*E(max(j-Nwindow+1,1):j)';  
      
        if(j >= Nwindow)
            
            
            if(abs(E(j))>reconstruction_error)
                
                YY =  Y(j-Nwindow+1:j);
                
                %startbacktime = k-Nwindow*interN;
                %xhat(:,startbacktime) = fmincon('myfun1_base1',xhat(:,startbacktime),[],[],[],[],lb,ub,[],options);
                %xhat(:,k) = fmincon('myfun1_base1',xhat(:,k)',[],[],[],[],lb,ub,[],options);
                actualk = k;
                xhat(:,k) = fminsearch('myfun1_base1',xhat(:,k));
                
                %Re-evaluate the estimated state (and all other
                %variables) from startbacktime to k with the new value
                %provided by fmincon
                jj = j-Nwindow;
%                for kk = startbacktime+1:k,
                   
%                    xhat(:,kk) = xhat(:,kk-1) + Ts*([-xhat(2,kk-1);xhat(1,kk-1)]  );%change
%                     yhat(kk) = xhat(1,kk); %change
%                    
%                     e(kk) =  y(kk)-yhat(kk);
%                     if(mod(kk,interN)==1)
%                         Y(jj) = y(kk);
%                         E(jj) = e(kk);
%                         J(jj) = E(max(j-Nwindow+1,1):j)*E(max(j-Nwindow+1,1):j)';  
%                         jj = jj + 1;
%                     end    
%               end
                
            end
         
        end
        j = j + 1;
    end
    
   
end
          
%%
figure(1)
ax(1) = subplot(4,1,1);
%plot(time,y,time,yhat,time,d,time,dhat,'LineWidth',2);
plot(time,y,time,yhat,'LineWidth',2);
grid on
%legend('y','\hat y','d','\hat d','FontSize',myfontsize);
ax(2) = subplot(4,1,2);
plot(time,e,'LineWidth',2)
grid on
ylabel('$e$','Interpreter','Latex','FontSize',myfontsize);
ax(3) = subplot(4,1,3);
plot(time(1:interN:end),J,'LineWidth',2);
grid on
ylabel('$J$','Interpreter','Latex','FontSize',myfontsize);
ax(4) = subplot(4,1,4);
plot(time,x,'LineWidth',2)
grid on
hold on
plot(time,xhat,'--','LineWidth',2)
ylabel('x','Interpreter','Latex','FontSize',myfontsize);
linkaxes(ax,'x')

%%
figure(2)
semilogy(time,abs(e))
grid on
ylabel('$\log (|e|)$','Interpreter','Latex','FontSize',myfontsize);
%%
figure(3)
plot(time,xhat-x)
grid on
ylabel('error','FontSize',myfontsize);
toc