function [a,siga,m,C,rms,chi2] = fitTS(t,e,n,sig_e,sig_n,toff,teq)
%
% function [a,siga,m,C,rms,chi2] = fitTS(t,e,n,sig_e,sig_n,toff,teq)
%
% Fits time series of horizontal GPS displacements for secular and
% periodic trends, coseismic offsets, and postseismic transients using
% nonlinear least squares.
%
% Included in the secular trends is an acceleration (quadratic) term.
%
% Input:
%   t   = vector of epochs (in decimal years)
%   e   = vector of displacements in east component (in mm)
%   n   = vector of displacements in north component (in mm)
%   sig_e = nominal displacement error in east component (in mm)
%   sig_n = nominal displacement error in north component (in mm)
%   toff = epochs of non-earthquake offsets (e.g. antenna changes)
%   teq = epochs of coseismic offsets (eartquakes)
%
% Output:
%   a   = estimated horizontal acceleration
%       a.ae (a.an) = value of acceleration in east (north) component
%       a.sigmae (a.sigman) = standard deviation in east (north)
%   m   = model vectors with all estimated parameters
%   C   = model covariance matrices for east and north
%   rms = root-mean-squared errors for eart and north
%   chi2 = reduced chi-squared values for east and north
%
%
% Andreas Mavrommatis, 2013

% --------------------------- Initialize -----------------------------

Noff = length(toff); Neq = length(teq);

% Initial guess
if 0 % if nonlinear parameters are free
    m0 = ones(7+Noff+3*Neq,1);
elseif 0 % if both nonlinear parameters (controling postseismic time scale) are fixed
    m0 = ones(7+Noff+Neq,1);
else % if only one nonlinear parameter is fixed
    m0 = ones(7+Noff+2*Neq,1);
    for k = 1:Neq
        m0(7+Noff+2*k) = 40;  % initial guess for nonlinear parameter
    end
end

% Lower and upper bounds on m - unconstrained
lb = -inf*ones(size(m0));
ub = inf*ones(size(m0));

% Bounds on nonlinear parameter
for k = 1:Neq
    lb(7+Noff+2*k) = 15;
    ub(7+Noff+2*k) = 60;
end

% ------------------------------ Solve --------------------------------

% Estimate model parameters, me and mn, using lsqcurvefit
% Levenberg-Marquardt and Gauss-Newton algorithms do not handle bound
% constraints; use default trust-region-reflective algorithm instead.

% Solve nonlinear unweighted least squares problem
[me,~,~,~,~,~,Je]  = lsqcurvefit(@(m,t) tsmodel(m,t,toff,teq), ...
    m0, t, e, lb, ub);
[mn,~,~,~,~,~,Jn]  = lsqcurvefit(@(m,t) tsmodel(m,t,toff,teq), ...
    m0, t, n, lb, ub);

% Compute model covariance matrices from Jacobians at solution
Ce = sig_e^2*((Je'*Je) \ eye(length(me)));
Cn = sig_n^2*((Jn'*Jn) \ eye(length(mn)));

% Compute residuals and weighted residuals
ehat = tsmodel(me,t,toff,teq);
nhat = tsmodel(mn,t,toff,teq);
rese = e - ehat;               resn = n - nhat;
reswe = rese/sig_e; reswn = resn/sig_n;

% ------------------------- Populate output ---------------------------

m.me = me;
m.mn = mn;

C.Ce = Ce;
C.Cn = Cn;

rms.rmse = sqrt(rese'*rese/length(t));
rms.rmsn = sqrt(resn'*resn/length(t));

% Reduced chi squared
chi2.chi2e = reswe'*reswe/(length(t) - length(me));
chi2.chi2n = reswn'*reswn/(length(t) - length(mn));

a.ae = me(3);
a.an = mn(3);

siga.sigae = sqrt(Ce(3,3));
siga.sigan = sqrt(Cn(3,3));

end