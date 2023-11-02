function ...
    [phase_support,...
    von_mises_pdf]=...
    make_von_mises_pdf(...
    mu,...
    kappa,...
    number_of_phase_bins)

% function ...
%     [phase_support,...
%     von_mises_pdf]=...
%     make_von_mises_pdf(...
%     mu,...
%     kappa,...
%     number_of_phase_bins);
%
% inputs --
%   mu: real-valued number in [-pi,pi), center of distribution
%   kappa: real-valued number in [0,inf), concentration of distriubtion
%       kappa=0 is flat (uniform) distribution on [-pi,pi)
%       kappa>0 is concentrated around mu, with larger kappa indicating
%       tighter concentration
%   number_of_phase_bins: integer number of values at which to evaluate the cdf
%
% outputs -- 
%   phase_support: 1 x number_of_phase_bins vector 
%                   evenly covers [-pi,pi)
%   von_mises_pdf: probability density function 
%           evaluated at the points in phase_support
%
% based on:
% http://en.wikipedia.org/wiki/Von_Mises_distribution

phase_support = (-pi:...
    2*pi/number_of_phase_bins:...
    pi-2*pi/number_of_phase_bins)';

von_mises_pdf = (2*pi*besseli(0,kappa))^-1*...
    exp(kappa*cos(phase_support-mu));
