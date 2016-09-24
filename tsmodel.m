
function F = tsmodel(m, t, toff, teq)
% Computes modeled time series data given vector of model parameters (m)
% and times (t). 

% Secular + periodic terms
F =  m(1) + m(2)*(t - t(1)) + 0.5*m(3)*(t-t(1)).^2 + ...
     m(4)*cos(2*pi*(t - t(1))) + m(5)*sin(2*pi*(t - t(1))) + ...
     m(6)*cos(4*pi*(t - t(1))) + m(7)*sin(4*pi*(t - t(1)));
 
% Add offsets 
N = length(toff);
for i = 1:N
    F = F + (t>=toff(i)).*m(7+i);
end

% Add transients for each earthquake - use solution for spring slider
% velocity-strengthening friction including loading velocity
for i = 1:length(teq)

if 0
    % full form (3 unknowns)
    B = m(7 + N + 3*i - 2); % amplitude
    a = m(7 + N + 3*i - 1); % sigma(a - b)/k
    b = m(7 + N + 3*i );    % v^inf
    
    F = F + B*(a*log(exp((b/a)*(t-teq(i))).*(t >= teq(i)) + exp(-1/a) - 1) ...
                - b*(t-teq(i)).*(t >= teq(i)) + 1).*(t >= teq(i));
   
else
    % reduced form (2 unknowns)
    B = m(7 + N + 2*i - 1); % amplitude
    a = m(7 + N + 2*i); % time decay parameter
    b = 60;
    
    F = F + B*(a*log(exp((b/a)*(t-teq(i)).*(t >= teq(i))) + exp(-1/a).*(t >= teq(i)) - 1.*(t >= teq(i))) ...
                - b*(t-teq(i)).*(t >= teq(i)) + 1.*(t >= teq(i))).*(t >= teq(i));
end 
            
end



