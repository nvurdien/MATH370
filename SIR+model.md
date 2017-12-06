



<style>
table {float:left}
</style>



# SIR Model Parameter list

 Parameter | What is it?                          | Baseline inputs| High V and Low M inputs | Low V and High M inputs
:----------|:-------------------------------------|:--------------:|:--------------:|:-------------:
 S         | Susceptive Population                | 1056           | 1056           | 1056
 I         | Infected Population                  | 322            | 322            | 322
 R         | Recovered Population                 | 0              | 0              | 0
 beta      | β -  infection rate                  | 0.00058        | 0.00058        | 0.00058
 gamma     | γ - recovery rate                    | 0.5            | 0.5            | 0.5
 v         | vaccination rate                     | 0.0            | 0.7            | 0.3
 m         | medication rate                      | 0.0            | 0.3            | 0.7
 t         | time                                 | 0 to 60        | 0 to 60        | 0 to 60
 y         | SIR rates in a y vector              | N/A            | N/A | N/A
 zeroloc   | checks when infection rate (I) is 0  | N/A            | N/A | N/A
 

## General Ebola SIR Model Function


```matlab
%%file SIRwithODE45.m
%SIR with medicine and vaccines

function [t, y] = SIRwithODE45(S, I, R, beta, gamma, v, m, fig)

y0 = [S I R];
tspan = [0 60];

% S is y(1), I is y(2)
[t, y] = ode45(@(t,y) [-(beta*y(1)*y(2)) - v*y(1); 
    beta*y(1)*y(2) - gamma*y(2) - m*y(2); 
    gamma*y(2) + v*y(1) + m*y(2)], tspan, y0);
% Plot the solutions for y against t.
figure(fig)
plot(t,y(:,1),'-',t,y(:,2),'-',t,y(:,3),'-')
title(sprintf('v = %.2f  m = %.2f', v, m));
xlabel('Time (Days)');
ylabel('Subpopulation Size (%)');

% Find when infected population becomes 0
zeroLoc = 1;
for i = 1:size(y,1)
    if y(i,2) < 0.5
        zeroLoc = i;
        break;
    end
end

hold on
plot(t(zeroLoc),y(zeroLoc,2), 'o')
legend('Susceptible','Infected','Recovered', '0.01% Infected')
hold off
```
