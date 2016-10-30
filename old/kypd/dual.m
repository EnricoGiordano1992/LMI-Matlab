function [X,Z,t,flag]=dual(d,H,G,n_basis,ntot_basis,matrix_info,options,F0,G0)
  
bigSize = 50000;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%30;
  N=matrix_info.N;
K=matrix_info.K;
nm=matrix_info.nm;

z = sdpvar(ntot_basis,1);



t = 0;
F = set([]);
for i=1:N
    Z{i}=F0{i};
    used_variables = (1:n_basis(i))+t;
    Hvec = reshape([H{used_variables}],nm(i)^2,[]);
    Z{i} = Z{i} + reshape(Hvec*z(used_variables),nm(i),nm(i));
    t=t+n_basis(i);
    F = F+set(Z{i}>0);
end;

F=F+set(G'*z+G0'==matrix_info.c);

if options.kypd.lowrank == 2
  options.sdpt3.smallblkdim=bigSize;
  foo = lmiinfo(F);
  base = size(foo.lin, 1) +...
	 size(foo.equ, 1) +...
	 size(foo.soc, 1) +...
	 size(foo.rlc, 1);
  small_ind = find(matrix_info.n<bigSize);
  big_ind = find(matrix_info.n>=bigSize);  
  foo = cell(length(n_basis), 1);
  tmp = 0; 
  for i=1:length(n_basis) %number of constraints
      tmp = tmp(end)+1:tmp(end)+n_basis(i);
     
        block_end = matrix_info.block_info(i);
        n = matrix_info.n(i);
        nm = matrix_info.nm(i);
        foo(i) = derankcross_sort(H(tmp), n_basis(i),block_end,n,nm);
  end
  
  if ~isempty(big_ind)
    for i = base+1:base+length(big_ind)
      options.sdpt3.schurfun{i}=@schurcompose;
    end
    options.sdpt3.schurfun_par(base+(1:length(big_ind)), 1) = foo(big_ind); 
    base = base + length(big_ind);
  end
  if ~isempty(small_ind)
    options.sdpt3.schurfun{base+1}=@schurcompose;
    options.sdpt3.schurfun_par{base+1, 1} = vertcat( foo{small_ind} );
    base = base + 1;
  end

end % if 'block'

options.solver = options.kypd.solver;

if options.kypd.reduce==1
   options.removeequalitites = 1;
   sol = solvesdp(F,ip(d,z),options);
else
   options.removeequalitites = 0;    
   sol = solvesdp(F,ip(d,z),options);
end;

flag = 0;
if sol.problem
    disp('KYPD terminated abnormally');
    disp('Error message received from YALMIP:');
    disp(sol.info);
    flag = 1;
end


for i=1:matrix_info.N
  Z{i}=double(Z{i});
  X{i}=double(dual(F(i)));
end;

t = sol.solvertime;