function J = myfun1_base1(x)

global YY Nwindow actualk  Ts interN

Yhat = zeros(1,Nwindow);
xp = x
j = Nwindow-1;

Yhat(Nwindow) = xp(1); %change

%iterate the new state 
for k = actualk-1:-1:actualk-(interN*(Nwindow-1))
    %xp(k) = xp(k+1) - T*dxp(k);
    xp = xp - Ts*([-xp(2);xp(1)]  );%change, use u(k-1) if necessary
    
    if(mod(k,interN)==1)
        Yhat(j) = xp(1);%change
        j = j - 1;
    end
end

J = (YY-Yhat)*(YY-Yhat)' 