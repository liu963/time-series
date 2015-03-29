function [TSf, CME, m] = cmefilt(TS,TS2,type,nw,verbose)
% [TSf, CME, m] = cmefilt(TS,TS2,nw)
%
% Remove common-mode errors (CME) in GPS time series.
%
% Input:
%   TS, TS2 are cell arrays with the following fields:
%        TS{i}.t        = time vector of ith station
%        TS{i}.site     = site id (string) of ith station
%        TS{i}.apcoords = a priori XYZ coordinates (3 x 1 vector)
%        TS{i}.X, TS{i}.Y, TS{i}.Z = X, Y, and Z coordinate time series
%    * TS should contain time series from *fiducial stations* (stations used to
%        estimate Helmert parameters for each epoch).
%    * TS2 should contain time series from *stations of interest* (stations upon
%        which the CME filter will be applied).
%   type = 'median' for median filter detrending, or
%          'lowpass' for low-pass filter detrending
%   nw   = window length (in days) if doing median filter, or
%          normalized cutoff frequency (from 0 to 1) if doing low-pass
%          filter
%
% Output:
%    TSf = cell array of filtered (CME removed) time series of stations of
%          interest, with the same structure as TS2.
%    CME = cell array of CME time series alone of stations of interest.
%          Same structure as TS2.
%    m   = (7 x Nt) vector of global Helmert parameters for each epoch (3
%          translations, 3 rotations, 1 scale), where Nt = number of unique
%          epochs
%
% Andreas Mavrommatis, 2013.

if nargin < 5
    verbose = 1;
end

% ------------------------------------------------------------------------
% STEP 1: Detrend raw XYZ time series of fiducial stations
% ------------------------------------------------------------------------

if verbose
    if strcmp(type,'median')
        disp('Running median filters...');
    elseif strcmp(type,'lowpass')
        disp('Running low-pass filters...');
    end
end

Ns = length(TS); % Number of fiducial stations

% Append antisymmetric images of original time series on boundaries in
% order to reduce boundary effects of the median filter
for i = 1:Ns;
    Nt = length(TS{i}.X);
    
    % Expand time vector
    TS{i}.t2(1:Nt) = TS{i}.t - (TS{i}.t(end) - TS{i}.t(1));
    TS{i}.t2(Nt+1:2*Nt) = TS{i}.t;
    TS{i}.t2(2*Nt+1:3*Nt) = TS{i}.t + (TS{i}.t(end) - TS{i}.t(1));
    
    % Antisymmetric images of X component
    TS{i}.Xl = -flipud(TS{i}.X) + 2*mean(TS{i}.X(1:15));
    TS{i}.Xr = -flipud(TS{i}.X) + 2*mean(TS{i}.X(end-15:end));
    
    % Append to time series (X component)
    TS{i}.X2(1:Nt) = TS{i}.Xl;
    TS{i}.X2(Nt+1:2*Nt) = TS{i}.X;
    TS{i}.X2(2*Nt+1:3*Nt) = TS{i}.Xr;
    
    % Now do the same for Y and Z
    TS{i}.Yl = -flipud(TS{i}.Y) + 2*mean(TS{i}.Y(1:15));
    TS{i}.Yr = -flipud(TS{i}.Y) + 2*mean(TS{i}.Y(end-15:end));
    
    TS{i}.Y2(1:Nt) = TS{i}.Yl;
    TS{i}.Y2(Nt+1:2*Nt) = TS{i}.Y;
    TS{i}.Y2(2*Nt+1:3*Nt) = TS{i}.Yr;
    
    TS{i}.Zl = -flipud(TS{i}.Z) + 2*mean(TS{i}.Z(1:15));
    TS{i}.Zr = -flipud(TS{i}.Z) + 2*mean(TS{i}.Z(end-15:end));
    
    TS{i}.Z2(1:Nt) = TS{i}.Zl;
    TS{i}.Z2(Nt+1:2*Nt) = TS{i}.Z;
    TS{i}.Z2(2*Nt+1:3*Nt) = TS{i}.Zr;
end

% Run and subtract median or low-pass filter from time series of fiducial sites

if strcmp(type,'lowpass'),
    Wn = nw;
    %     Wn = 0.001;
end
for i = 1:Ns
    Nt = length(TS{i}.X);
    
    if strcmp(type,'median'),
        TS{i}.X2f = medfilt1(TS{i}.X2,nw); % Median filter of extended time series
    elseif strcmp(type,'lowpass')
        TS{i}.X2f = lpfilter(TS{i}.X2,Wn,0); % Low pass of extended time series
    end
    TS{i}.Xf = TS{i}.X2f(Nt+1:2*Nt); % Extract median filter for original ts
    TS{i}.Xd = TS{i}.X' - TS{i}.Xf; % Detrend by subtracting median filter
    
    if strcmp(type,'median'),
        TS{i}.Y2f = medfilt1(TS{i}.Y2,nw);
    elseif strcmp(type,'lowpass')
        TS{i}.Y2f = lpfilter(TS{i}.Y2,Wn,0); % Low pass of extended time series
    end
    TS{i}.Yf = TS{i}.Y2f(Nt+1:2*Nt);
    TS{i}.Yd = TS{i}.Y' - TS{i}.Yf;
    
    if strcmp(type,'median'),
        TS{i}.Z2f = medfilt1(TS{i}.Z2,nw);
    elseif strcmp(type,'lowpass')
        TS{i}.Z2f = lpfilter(TS{i}.Z2,Wn,0); % Low pass of extended time series
    end
    TS{i}.Zf = TS{i}.Z2f(Nt+1:2*Nt);
    TS{i}.Zd = TS{i}.Z' - TS{i}.Zf;
end



% ------------------------------------------------------------------------
% STEP 2: Estimate global daily Helmert parameters from fiducial stations
% ------------------------------------------------------------------------

if verbose, disp('Estimating daily Helmert parameters...'); end

% Make Helmert transformation matrices for each fiducial station
H = cell(Ns,1);
for i = 1:Ns
    H{i} = make_helmert([TS{i}.apcoords]);
end

% Find unique epochs from the set of fiducial stations
epochs = cell(Ns,1);
for k = 1:Ns
    epochs{k} = TS{k}.t;
end
all_epochs = cat(1,epochs{:});
unique_epochs = unique(all_epochs);

% Do for every unique epoch
m = zeros(7,length(unique_epochs)); % Helmert parameters for each unique epoch
for k = 1:length(unique_epochs)
    ti = unique_epochs(k); % current epoch
    d = []; G = [];
    % Do for every fiducial station
    for i = 1:Ns
        j = find(TS{i}.t == ti);  % epoch index in the current station
        if ~isempty(j)   % if this epoch does not exist, don't include this station
            d = [d; TS{i}.Xd(j); TS{i}.Yd(j); TS{i}.Zd(j)]; % add to data vector
            G = [G; H{i}]; % add to stacked Helmert matrix
            m(:,k) = G\d; % estimate Helmert parameters for this day
        end
    end
end

if verbose, disp('Helmert parameter estimation done.'); end


% --------------------------------------------------------------------------
% STEP 3: Apply Helmert transformation on each epoch to stations of interest
% --------------------------------------------------------------------------

% For each station of interest, see if the station contains the current
% unique epoch. If yes, compute predicted common-mode position for that
% epoch, j, d_cme(j) = H*m(k), where H = Helmert matrix of this station,
% and m(k) = Helmert parameters for absolute epoch k.

Ns2 = length(TS2); % Number of stations of interest

% Make Helmert transformation matrices for each station of interest
H2 = cell(Ns2,1);
for i = 1:Ns2
    H2{i} = make_helmert([TS2{i}.apcoords]);
end

% Compute predicted CME
CME = cell(Ns2,1);

% Do for every station of interest
for i = 1:Ns2
    
    CME{i}.t = TS2{i}.t;
    CME{i}.site = TS2{i}.site;
    CME{i}.apcoords = TS2{i}.apcoords; % This fields are included only for reference
    
    % Do for every unique epoch
    for k = 1:length(unique_epochs)
        ti = unique_epochs(k); % current epoch
        j = find(TS2{i}.t == ti);  % epoch index in the current station
        if ~isempty(j)   % if this epoch does not exist, skip this station
            dcm(:,j) = H2{i}*m(:,k); % update predicted CME for this epoch
            CME{i}.X = dcm(1,:)'; % Update CME vectors each time
            CME{i}.Y = dcm(2,:)';
            CME{i}.Z = dcm(3,:)';
        end
    end
    clear dcm % Restart for next station
    
    if verbose,
        disp(['CME filter applied on station ',num2str(TS2{i}.site),', ',...
            num2str(i),' out of ',num2str(Ns2)]);
    end
end

if verbose,
    disp('Filtering...');
end

% Substract predicted CME to get filtered time series, TSf
TSf = cell(Ns2,1);
for i = 1:Ns2
    
    TSf{i}.t = TS2{i}.t;
    TSf{i}.site = TS2{i}.site;
    TSf{i}.apcoords = TS2{i}.apcoords;
    
    TSf{i}.X = TS2{i}.X - CME{i}.X;
    TSf{i}.Y = TS2{i}.Y - CME{i}.Y;
    TSf{i}.Z = TS2{i}.Z - CME{i}.Z;
    
end


if verbose, disp('All done!'); end

end

% ------------------------------------------------------------------------

function H = make_helmert(XYZ)
% H = make_helmert(XYZ)
%
% Make Helmert matrix H, so that d = Hm, where d = observed position and m
% = unknown model parameters = [s c1 c2 c3 r1 r2 r3], where s = scale, ci =
% translations, and ri = rotations.
%
% Input: XYZ = 3*Nsta x 1 vector of (geocentric) a priori coordinates of
%              selected stations
%
% APM, 2013

Nsta = length(XYZ)/3;               % Number of stations

% Rotation matrix
R = zeros(Nsta,3);
for i = 1:Nsta
    X = XYZ(i); Y = XYZ(i+1); Z = XYZ(i+2); % Position components
    R(3*i-2:3*i,:) = [0 Z -Y; -Z 0 X; Y -X 0]; % Rotation matrix
end

H = [XYZ repmat(eye(3),Nsta,1) R];
% scale    translation     rotation

end

% ------------------------------------------------------------------------

function xf = lpfilter(x,Wn)
% xf = lpfilter2(x,Fs,Fc)
% 
% Low-pass filter.
% Inputs: x = data vector
%         Wn = normalized cutoff frequency (between 0 and 1, where 1
%         corresponds to the Nyquist frequency)

[B,A] = butter(1,Wn);   % Butterworth filter
xf = filtfilt(B,A,x);

end
