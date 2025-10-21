% Define the model function 
function J_model = Zhang_model(p)
    global V Voc Jsc T k q J;
     n = p(1);
     Rs = p(2);
     Rsh = p(3);
     J_cal = ((n*k*T)/(q*Rs))*lambertw((q*Rs)/(n*k*T)*(Jsc-Voc/(Rs+Rsh))*(exp((-q*Voc)/(n*k*T)))*(exp(q/(n*k*T)*(Rs*Jsc+((Rsh.*V)./(Rsh+Rs))))))+V./Rs-Jsc-((Rsh.*V)./(Rs*(Rsh+Rs)));
    J_model = J_cal-J;
end
