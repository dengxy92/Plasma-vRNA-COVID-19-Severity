%% Simulation of SIV model wtih delay as a DDE, and logistic target cell growth
function [sol,p] = simulation_virus_model_with_delay_no_tinf(p,tspan)

%Run integration again until completion
sol = ode15s(@odefun,tspan,p.IC);

%------------------------------------------------------------------------
function dydt = odefun(t,y)

T = y(1);
I = y(2);
V = y(3);

dT = -p.bet*T*V;%p.r*T*(1-(T+I1+I2)/p.T_max)
dI = p.bet*T*V-p.d_I*I;
dV = p.p*I-p.d_V*V; %infectious
 
dydt = [dT;dI;dV];

end

%-------------------------------------------------------------------------
end
