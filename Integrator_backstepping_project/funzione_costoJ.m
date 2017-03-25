function J = funzione_costoJ(x)

global YY Nwindow actualk  Ts interN

Yhat = zeros(1,Nwindow);
xp = x;
j = Nwindow-1;

Yhat(Nwindow) = xp(1); %change

%iterate the new state 
for k = actualk-1:-1:actualk-(interN*(Nwindow-1))
   
    xp(1) = xp(1) + Ts*xp(2);
    xp(2) = xp(2) + Ts*(-xp(3));
    xp(3) = xp(3) + Ts*(-0.7*xp(3) - 1*xp(2) + xp(1) - xp(1)^3);
    
    if(mod(k,interN)==1)
        Yhat(j) = xp(1);%change
        j = j - 1;
    end
end

J = (YY-Yhat)*(YY-Yhat)' ;