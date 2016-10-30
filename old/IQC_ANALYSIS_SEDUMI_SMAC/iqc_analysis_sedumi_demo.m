%
% Demo of iqc_analysis_sedumi.m on an LFT transport aircraft model
%
clear;
%
% *** Closed loop LFT model of the transport aircraft ***
%
% I/O of sys_M:
% 1   : deadzone nonlinearity, corresponding to a rate limiter
% 2:25: repeated parameters CG, Mach, Mass, Conventional Airspeed Vc (see blk for the degree of repetition)
% 26: performance block, corresponding to a transfer function between the reference input qc and the tracking error
%
load iqc_analysis_sedumi_demo_data.mat; % sys_M and blk are loaded
%
% CG and Mass are considered as LTI parameters
% Mach and Vc are LTI or LTV parameters, with or without a bounded rate of variation
%
fprintf('\n');
while 1
  type_param=input(' LTV (0), LTI (1) or bounded rate (2) Mach and Vc parameters (default value is 0): ');
  if isempty(type_param) type_param=0; end
  if any(type_param==[0 1 2]) break; end
end
fprintf('\n');
%
delta=[];
delta(1).dim=1;
delta(1).type='nl-popov';
delta(1).kmax=0.3;
%
for i=1:size(blk,1)
  delta(i+1).dim=-blk(i,1);
  if any(i==[1 3])
    delta(i+1).type='param-lti';
    delta(i+1).poles=2;
  elseif type_param==0
    delta(i+1).type='param-ltv';
  elseif type_param==1
    delta(i+1).type='param-lti';
    delta(i+1).poles=2;        
  else
    delta(i+1).type='param-bounded-rate';
    delta(i+1).poles=2;
    delta(i+1).rate_max=1;
  end
end
%
options=[];
options.visu=1;
options.puls=linspace(0,10,100);
options.method=1; % 2 if using the KYPD solver
tic;
[gain,info]=iqc_analysis_sedumi(sys_M,delta,options);
fprintf('\n');
toc;
%
fprintf('\n Robust performance level = %6.3f\n',gain);
fprintf('\n Info provided by sedumi.m and iqc_analysis_sedumi.m:\n');
disp(info);

