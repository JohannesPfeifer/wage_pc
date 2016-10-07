%/*************************************************************************
% * Replication file for Born/Pfeifer (2016): "The New Keynesian Wage
% * Phillips Curve: Calvo vs. Rotemberg"
% * 
% * Computes the Rotemberg wage adjustment cost parameter corresponding
% * to a particular average wage duration
% *
% ************************************************************************\
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Use this to compute the total elasticity for multiplicatively separable preferences
%
% sigma=2
% eta=0.34
% labor_share_in_leisure=1/3/(2/3);
% eps_tot=(1-((1-eta)*(sigma-1))/(eta*(1-sigma)-1))*labor_share_in_leisure

EHL_dummy = 1;              % set to 1 if EHL-framework, 0 for SGU
fixed_cost_dummy = 1;       % set to 1 if output to which Xi refers is net of fixed costs from monopolistic pricing power
par.MRS_elasticity = 1;     % required if EHL framework, can be set to any value otherwise
par.epsilon_p = 11;         % elasticity of substitution in goods market, determines markup;  required if fixed costs are to be subtracted, can be set to any value otherwise 
par.epsilon_w = 11;         % elasticity of substitution labor services
par.theta_w = 0.75;         % target price duration
par.betta = 0.99;           % discount factor
par.tau_n = 0.21;           % steady state labor tax
par.alpha = 1/3;            % capital in production function

[phi_Rotemberg, Rotemberg_slope]=compute_Rotemberg_parameter(par,fixed_cost_dummy,EHL_dummy);


fprintf('\nSlope: %8.6f\n',Rotemberg_slope)
fprintf('phi: %8.6f\n',phi_Rotemberg)

%% Make tables
labels=strvcat('SGU','EHL');
par_old=par;
column_iter=1;
headers=' ';
headers_tex=' ';
epsilon_w_vec=[6 11 21];
for par_iter=1:length(epsilon_w_vec)
    for row_iter=1:2
        EHL_dummy=row_iter-1;
        par.epsilon_w=epsilon_w_vec(par_iter);
        [phi_Rotemberg(row_iter,column_iter)]=compute_Rotemberg_parameter(par,fixed_cost_dummy,EHL_dummy);
    end
    headers=strvcat(headers,['eps_w=',num2str(epsilon_w_vec(par_iter))]);
    headers_tex=strvcat(headers_tex,['\varepsilon_w=',num2str(epsilon_w_vec(par_iter))]);
    column_iter=column_iter+1;
end
par=par_old; %reset baseline
lh = size(labels,2)+2;
if ~exist('Tables','dir')
    mkdir('Tables')
end
make_latex_table('Par_table','','substitution_elasticity',headers_tex,labels,phi_Rotemberg(:,end-length(epsilon_w_vec)+1:end),lh,4,2);
clear epsilon_w_vec EHL_dummy

tau_n_vec=[0 0.21 0.4];
for par_iter=1:length(tau_n_vec)
    for row_iter=1:2
        EHL_dummy=row_iter-1;
        par.tau_n=tau_n_vec(par_iter);
        [phi_Rotemberg(row_iter,column_iter)]=compute_Rotemberg_parameter(par,fixed_cost_dummy,EHL_dummy);
    end
    headers=strvcat(headers,['tau_n=',num2str(tau_n_vec(par_iter))]);
    headers_tex=strvcat(headers_tex,['\tau^n=',num2str(tau_n_vec(par_iter))]);
    column_iter=column_iter+1;
end
par=par_old; %reset baseline
make_latex_table('Par_table','','labor_tax',headers_tex([1,end-length(tau_n_vec)+1:end],:),labels,phi_Rotemberg(:,end-length(tau_n_vec)+1:end),lh,4,2);
clear tau_n_vec EHL_dummy

make_latex_table('Par_table','','first_part',headers_tex,labels,phi_Rotemberg,lh,4,2);

first_part_table_end=size(headers_tex,1);


MRS_n_vec=[0.25 1 1.5];
for par_iter=1:length(MRS_n_vec)
    for row_iter=1:2
        EHL_dummy=row_iter-1;
        par.MRS_elasticity=MRS_n_vec(par_iter);
        [phi_Rotemberg(row_iter,column_iter)]=compute_Rotemberg_parameter(par,fixed_cost_dummy,EHL_dummy);
    end
    headers=strvcat(headers,['epsilon_n^MRS=',num2str(MRS_n_vec(par_iter))]);
    headers_tex=strvcat(headers_tex,['\varepsilon_n^{MRS}=',num2str(MRS_n_vec(par_iter))]);
    column_iter=column_iter+1;
end
par=par_old; %reset baseline
make_latex_table('Par_table','','MRS_elasticity',headers_tex([1,end-length(MRS_n_vec)+1:end],:),labels,phi_Rotemberg(:,end-length(MRS_n_vec)+1:end),lh,4,2);
clear MRS_n_vec EHL_dummy

betta_vec_vec=[0.985 0.99 0.995];
for par_iter=1:length(betta_vec_vec)
    for row_iter=1:2
        EHL_dummy=row_iter-1;
        par.betta=betta_vec_vec(par_iter);
        [phi_Rotemberg(row_iter,column_iter)]=compute_Rotemberg_parameter(par,fixed_cost_dummy,EHL_dummy);
    end
    headers=strvcat(headers,['beta=',num2str(betta_vec_vec(par_iter))]);
    headers_tex=strvcat(headers_tex,['\beta=',num2str(betta_vec_vec(par_iter))]);
    column_iter=column_iter+1;
end
par=par_old; %reset baseline
make_latex_table('Par_table','','discount_factor',headers_tex([1,end-length(betta_vec_vec)+1:end],:),labels,phi_Rotemberg(:,end-length(betta_vec_vec)+1:end),lh,4,2);
clear betta_vec_vec EHL_dummy

make_latex_table('Par_table','','second_part',headers_tex([1,first_part_table_end+1:end],:),labels,phi_Rotemberg(:,first_part_table_end:end),lh,4,2);

fixed_cost_dummy=0;

for row_iter=1:2
    EHL_dummy=row_iter-1;
    [phi_Rotemberg(row_iter,column_iter)]=compute_Rotemberg_parameter(par,fixed_cost_dummy,EHL_dummy);
    [phi_Rotemberg_fixed(row_iter,1)]=compute_Rotemberg_parameter(par,1,EHL_dummy);
    [phi_Rotemberg_fixed(row_iter,2)]=compute_Rotemberg_parameter(par,0,EHL_dummy);
end

make_latex_table('Par_table','','fixed_cost',strvcat(' ','Baseline','No fixed cost'),labels,phi_Rotemberg_fixed,lh,4,2);
