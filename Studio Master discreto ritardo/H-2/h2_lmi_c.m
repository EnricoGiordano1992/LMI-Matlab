function out = h2_lmi_c(A,B,C,param)
% function out = h2_lmi_c(A,B,C,param)
%
% Evaluate the H-2 norm of a continuous-time linear system. The LMIs are programmed using YALMIP and can be 
% solved by any LMI solver supported by YALMIP (SeDuMi is the default). 
%
% input:  (A,B,C) -> state-space matrices
%        
% output: out.h2       -> H-2 norm (0 if unfeasible)
%         out.cpusec   -> cpu time to solve the LMIs (seconds)
%         out.cpusec_m -> cpu time to mount the LMIs (seconds)
%         out.P        -> Lyapunov matrix
%	      out.V        -> number of decision variables
%	      out.L        -> number of LMIs rows
%
% Date: 04/09/2009
% Author: ricfow@dt.fee.unicamp.br

if nargin == 4
    if isfield(param,'h2')
        h2 = param.h2;
    else
        h2 = 0;
    end
    if isfield(param,'precision')
        precision = param.precision;
    else
        precision = 1e-7;
    end
else
    h2 = 0;
    precision = 1e-7;
end

%determine the number of vertices of polytope
order = size(A,1);
inputs = size(B,2);
outputs = size(C,1);

out.cpusec_m = clock;

%new LMI system
LMIs = set([]);
%LMI rows counter
out.L = 0;

%create the variables
if h2 == 0
    mu = sdpvar(1);
    obj = mu;
else
    mu=h2*h2;
    obj = [];
end

P = sdpvar(order,order,'symmetric');
LMIs = [LMIs, P > 0];
out.L = out.L + order;
    
%trace condition (mu > trace( B' P B) )
LMIs = [LMIs, mu > trace(B'*P*B)];

%gramian condition (A P + P A' + C'C < 0)
LMIs = [LMIs, A'*P + P*A + C'*C < 0];

out.L = out.L + order + 1;

out.cpusec_m = etime(clock,out.cpusec_m);

out.V = size(getvariables(LMIs),2);

sol = solvesdp(LMIs,obj,sdpsettings('verbose',0,'solver','sedumi'));
out.cpusec = sol.solvertime;
p=min(checkset(LMIs));

out.h2 = 0;
%capturing the solutions (if ones exist)
if p > -precision
    out.h2 = sqrt(double(mu));    
    out.P = double(P);
end