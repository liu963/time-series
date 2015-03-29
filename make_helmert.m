function [H,R] = make_helmert(XYZ)
% H = make_helmert(XYZ)
%
% Make Helmert matrix H, so that d = Hm, where d = observed position and m
% = unknown model parameters = [s c1 c2 c3 r1 r2 r3], where s = scale, ci =
% translations, and ri = rotations.
% 
% Input: XYZ = 3*Nsta x 1 vector of (geocentric) a priori coordinates of 
%              selected stations
%
% Outputs: H = Helmert matrix [scale translation rotation]
%          R = Rotation matrix (last 3 columns of H)
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
 
 
 
    
