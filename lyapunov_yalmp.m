%This example illustrates the definition and solution of a simple
%semidefinite programming problem. 

%Given a linear dynamic system x' = Ax, our goal is to prove stability by
%finding a symmetric matrix P satisfying 

%Define a stable matrix A and symmetric matrix P (remember: square matrices
%are symmetric by default) 
A = [-1 2 0;-3 -4 1;0 0 -2];
P = sdpvar(3,3);

%Having P, we are ready to define the constraints.
F = [P >= 0, A'*P+P*A <= 0];

%To avoid the zero solution or an unbounded solution, we constrain the
%trace of the matrix (Of course, this is not the only way. We could have
%used, e.g., the constraint P>I instead)  
F = [F, trace(P) == 1];

%At this point, we are ready to solve our problem. But first, we display
%the collection of constraints to see what we have defined. 
F

%We only need a feasible solution, so one argument is sufficient when we
%call optimize to solve the problem. 
optimize(F);
Pfeasible = value(P);

%The resulting constraint satisfaction is easily investigated with
%checkset. 
checkset(F)

%Minimizing, e.g., the top-left element of P is done by specifying an
%objective function. 
F = [P >= 0, A'*P+P*A <= 0, trace(P)==1];
optimize(F,P(1,1));

%We can easily add additional linear inequality constraints. If we want to
%add the constraint that all off-diagonal elements are larger than zero,
%one approach is (remember, standard MATLAB indexing applies)  
F = [P >= 0, A'*P+P*A <= 0, trace(P)==1, P([2 3 6])>=0];
optimize(F,P(1,1));

%Since the variable P([2 3 6]) is a vector, the constraint is interpreted
%as a standard linear inequality, according to the rules introduced in the
%basic tutorial.   