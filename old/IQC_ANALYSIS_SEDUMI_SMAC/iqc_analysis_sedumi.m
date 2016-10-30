function [gain,info]=iqc_analysis_sedumi(sys_M,Delta,options)
%
% THIS SOFTWARE IS ONLY INTENDED FOR STRICT ACADEMIC USE
%
% THE SEDUMI SOLVER (sedumi.m and associated routines) NEEDS TO BE AVAILABLE IN THE PATH
%
% IF USING THE KYPD SOLVER, YALMIP AND THE KYPD SOLVER (kypd_solver.m and associated routines) ALSO NEED TO BE AVAILABLE IN THE PATH
%
% call:  [gain,info]=iqc_analysis_sedumi(sys_M,Delta,options)
%
% Computation of a guaranteed value "gain" of the worst case norm of F_u(sys_M,Delta)
% sys_M = state-space description of M(s) (ss format)
% First I/O of sys_M associated to Delta, described by input argument "Delta"
%
% "info" = output of the Sedumi solver (type "help sedumi")
% To be checked in "info":
% 1/ Optimization problem considered as feasible if info.pinf = info.dinf = 0
% 2/ Serious numerical problems if info.numerr = 2
%
% Additional information in "info":
% info.nx: number of states of the augmented problem
% info.nbr_primal_variables: number of "primal" variables = number of variables of the primal state-space LMI problem
% info.nbr_dual_variables  : number of "dual"   variables
%
% options.perf = 1 --> Robust performance problem
%                      output argument "gain" = minimized guaranteed performance level
%                                             = Inf if the robust stability condition is not satisfied
%                2 --> Robust performance test
%                      Is the minimized guaranteed performance level less than options.L2gain ?
%                      output argument "gain" = options.L2gain or Inf
%                      Rk: gain = Inf may also mean that the robust stability test failed.
%                0 --> Robust stability problem (M(s) may contain extra I/O for performance)
%                      output argument "gain" = 0 or Inf, if the robust stability test passed or failed
% Default value for options.perf: 1
%
% If info.numerr = 2 is obtained for options.perf = 1, the minimized L2 gain is validated by performing a robust performance test,
% where the tested value corresponds to the minimized L2 gain relaxed by the tolerance options.epsL2gain, except if options.epsL2gain <= 0
% Default value for options.epsL2gain: 0.05 (5 percent)
%
% Only if options.perf = 1:
% info.reliability = 0 --> numerr = 0 at the output of the Sedumi solver (no numerical problem)
%                    1 --> numerr = 1 (slight numerical problems)
%                    2 --> numerr = 2 + relaxed L2 gain validated with a robust performance test
%                    3 --> Optimization problem considered as unfeasible
%                   -1 --> Robust L2 gain strictly less than the nominal one (obtained for Delta = 0)
%                   -2 --> numerr = 2 + relaxed L2 gain invalidated with a robust performance test
%
% Special case of options.perf = 3 --> Robust performance problem solved with robust performance tests and a dichotomy search
%                                      between options.L2gain_min and options.L2gain_max
% Default value for options.L2gain_min: L2 gain of the nominal model (corresponding to Delta = 0)
% options.L2gain_max should be specified.
%
% options.visu = 0/1 --> visualization of intermediate results of the Sedumi solver
% Default value: 0
%
% options.puls is a frequency gridding used to provide a "Bode" plot of the constraint for the optimal value of the multiplier.
% The less negative the constraint, the less satisfied !
% Default value: [] (no plot)
% Rk: visualizing the frequency-domain constraint can be useful for better choosing the multiplier poles.
%
% options.method = 1 --> Primal state-space problem solved with Sedumi
%                  2 --> Dual   state-space problem solved with Sedumi and the KYPD solver
% Default value: 1
% Rk1: options.method = 2 is available only for options.perf = 1
%      (it seems difficult to check the optimization problem feasibility with the KYPD solver)
% Rk2: output argument "info" = [] if options.method = 2
% Rk3: Yalmip is not necessary if options.method = 1, but necessary if options.method = 2
%
% *** DESCRIPTION OF DELTA ***
%
% Delta(i).dim = size of block # i of uncertainty.
%
% Delta(i).type='param-ltv' --> LTV uncertainty.
% Rk: Delta(i).magnitude may be specified as the magnitude of the uncertainty. Default value is 1.
%
% Delta(i).type='param-lti' --> LTI uncertainty, real poles of the multipliers in Delta(i).poles
% Rk1: positive values for the "poles" if stable multipliers.
% Rk2: Delta(i).magnitude may be specified as the magnitude of the uncertainty. Default value is 1.
%
% Delta(i).type='param-bounded-rate' --> LTV uncertainties with bounded rate Delta(i).rate_max,
%                                        real poles of the multipliers in Delta(i).poles
% Rk1: positive values for the "poles" if stable multipliers.
% Rk2: Delta(i).magnitude may be specified as the magnitude of the uncertainty. Default value is 1.
%
% Delta(i).type='nl-circle' --> SISO nonlinearity in a sector [0, Delta(i).kmax],  use of the circle criterion.
%
% Delta(i).type='nl-popov'  --> SISO nonlinearity in a sector [0, Delta(i).kmax], use of the Popov criterion.
%
% G. FERRERES, Onera/DCSD, 2015
%
if nargin~=3 error('IQC_ANALYSIS_SEDUMI: 3 input arguments'); end
%
[gain,info]=iqc_analysis_sedumi_code(sys_M,Delta,options);
