function output = hinf_norm_c_yal(A,B,C,D,param)
%function output = hinf_norm_c_yal(A,B,C,D,param)
%
% Compute the H-infininty norm of a continuous-time linear system using linear
% matrix inequalities. The LMIs are  programmed using YALMIP and can be solved 
% by any LMI solver supported by YALMIP (SeDuMi was used with the default options).
% inputs:  A,B,C,D   -> vertices of the polytope (cell array)
%          param   -> optional parameters
%          param.hinf -> if this parameter is different of zero, the routine will
%                        check the feasibility of the LMIs for this value of the hinf
%                        norm.
%          param.tol -> tolerance of the feasbility of the LMIs. 
%
% outputs: output.hinf     -> H-infinity guaranteed cost (0 if unfeasible)
%          output.cpusec   -> cpu time to solve the LMIs (seconds)
%          output.cpusec_m -> cpu time to mount the LMIs (seconds)
%          output.Pk       -> Solution variables P
%          output.K        -> number of scalar variables used in the optmization problem
%          output.L        -> number of LMI rows used in the optmization problem
%          output.delta    -> minimal primal residual returned by the LMI solver (SeDuMi is the default).
%
% Example: 2 states and 2 vertices
% A = [-0.9  0.2;
%       -0.5  -1.9 ];
% B = [1;0];
% C = [1 0];
% D = 0;
% output = hinf_norm_c_yal(A,B,C,D)
%
% Date: 29/08/2011
% Author: ricfow@dt.fee.unicamp.br

if nargin == 5
    if isfield(param,'hinf')
        hinf = param.hinf;
    else
        hinf = 0;
    end
    if isfield(param,'tol')
        tol = param.tol;
    else
        tol = 1e-7;
    end    
else
    hinf = 0;
    tol = 1e-7;
end

%determine the number of vertices of polytope
order = size(A,1);
inputs = size(B,2);
outputs = size(C,1);

output.cpusec_m = clock;

%new LMI system
LMIs = [];

if hinf == 0
    mu = sdpvar(1,1);
    obj = mu;
else
    mu = hinf^2;
    obj = [];
end

P = sdpvar(order,order,'symmetric');
LMIs = [LMIs, P > 0];

% Bounded real lemma  < 0
T11 = A'*P + P*A;
T12 = P*B;
T13 = C';
T22 = -eye(inputs);
T23 = D';
T33 = -mu*eye(outputs);
T = [T11 T12 T13;
    T12' T22 T23;
    T13' T23' T33];
LMIs = [LMIs, T < 0];

output.cpusec_m = etime(clock,output.cpusec_m);

% determine the number of LMI rows
output.L = 0;
for i=1:size(LMIs,1)
    output.L = output.L + size(LMIs{i},1);
end
% determine the number of scalar variables
output.V = size(getvariables(LMIs),2);

% solve the LMIs
sol = solvesdp(LMIs,obj,sdpsettings('verbose',0,'solver','sedumi'));
% evaluate the elapsed time to solve the LMIs set
output.cpusec_s = sol.solvertime;

% retrieving the minimal primal residual
p=min(checkset(LMIs));
output.delta = p;

output.feas = 0;
% capturing the solutions (if ones exist)
if p > -tol
    output.P = double(P);
    output.hinf = sqrt(double(mu));
    output.feas = 1;
end