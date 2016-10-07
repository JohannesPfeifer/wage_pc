function [phi_Rotemberg, Rotemberg_slope]=compute_Rotemberg_parameter(par,fixed_cost_dummy,EHL_dummy)
% function [phi_Rotemberg, Rotemberg_slope]=compute_Rotemberg_parameter(par,fixed_cost_dummy,EHL_dummy)
% Inputs
%   - par               [structure]     contains parameter values
%   - fixed_cost_dummy  [boolean]       indicator for fixed costs
%   - EHL_dummy         [boolean]       1: EHL, 0: SGU
%
% Outputs
%   - phi_Rotemberg     [scalar]        implied Rotemberg parameter 
%   - Rotemberg_slope   [scalar]        slope of the wage Phillips Curve

% Copyright (C) 2016 Johannes Pfeifer and Benjamin Born
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

%Compute cost share
if fixed_cost_dummy
    cost_share=(1-par.alpha);
else
    cost_share=(par.epsilon_p-1)/par.epsilon_p*(1-par.alpha);
end

Calvo_slope = ((1-par.theta_w)*(1-par.betta*par.theta_w))/(par.theta_w*(1+EHL_dummy*par.epsilon_w*par.MRS_elasticity));
phi_Rotemberg = (par.epsilon_w-1)*(1-par.tau_n)*cost_share*par.theta_w*(1+EHL_dummy*par.epsilon_w*par.MRS_elasticity)/((1-par.theta_w)*(1-par.betta*par.theta_w));
Rotemberg_slope = (par.epsilon_w-1)*(1-par.tau_n)*cost_share/phi_Rotemberg;