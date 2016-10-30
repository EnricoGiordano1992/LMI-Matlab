function [u,P,x,Z,soltime]=kypd_solver(matrix_info,options)
%
% [u,P,x,Z,soltime]=kypd_solver(matrix_info,options)
%
% Solves semidefinite programs originating from the KYP lemma. 
% If transform=1 or the system matrices are not Hurwitz feedback 
% will be applied on the system. If the number of variables can be
% reduced this will be done to increase efficiency. The solution
% produced is however a solution to the original formulation of
% the problem. 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    options = sdpsettings;
end

% Check that the Mij matrices are linearly independent, check 
% stabilizability and make the system Hurwitz. Reduce the system
% if necessary.

[matrix_info,temp,T,Pbar,V,c,L,checkM,checkstab]=...
                    check_and_transform(matrix_info,...
                    options.kypd.transform,options.kypd.rho,...
                    options.kypd.lyapunovsolver,options.kypd.lowrank);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve the problem

[u,P,x,Z,soltime,errorflag]=kypd(matrix_info,options);

% Transform the solutions (if block-diagonalization has been done)
if options.kypd.lowrank & ~errorflag
    for i = 1:matrix_info.N
        if matrix_info.n(i)
            P{i} = matrix_info.T_diag{i}.'\P{i}/matrix_info.T_diag{i};
            matrix_info.C{i} = ...
                matrix_info.T_diag{i}*matrix_info.C{i}*matrix_info.T_diag{i}.';
        end
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transform the solutions to satisfy the original LMIs.
if ~errorflag

    [P,x,Z]=transform_solution(temp,x,P,Pbar,Z,T,V,c,L,checkM,checkstab,...
        options.kypd.lyapunovsolver);

    c = matrix_info.c;
    tmp = 0;
    for i = 1:matrix_info.N
        if matrix_info.n(i) > 0
            for j = 1:matrix_info.K
                c(i) = c(i) + trace(matrix_info.C{i}*Pbar{i,j+1});
            end
            tmp = tmp + trace(matrix_info.C{i}*P{i});
        end
    end
    u = c.'*double(x)+tmp;
else
    fprintf('\n \n  The solver kypd_solver has terminated due to an error! \n \n');
end

