function M = schurcompose( W, Zinv, caux )
% M = schurcompose( W, Zinv, caux )
% 
% "schurfun" to be used with SDPT3.
%
% Input arguments:
%
% o  <W> is the scaling matrix.
%
% o  <Zinv> is ignored since it will be equal to W.
%
% o  <caux> are the auxillary parameters that we use to compute 
%    the Schur block.  It is a cell array of size [ 1 1 ], containing
%    a cell array of size [ <n_blk> 3 ], where <n_blk> is the number
%    of smaller "constraint blocks" in the block diagonal <M>.
%
% Output:
%
% o  <M> is the computed Schur block.
%
%
% The output is computed by computation and assembly of the individual
% blocks.  Each block is computed by the function lowrankschur6, which
% is called with the corresponding part of the scaling matrix and
% the corresponding row of <caux>{1} as arguments.
%
% See also: lowrankschur6

aux = caux{ 1 };

if size( aux, 1 ) == 1
  M = lowrankschur6( W, aux );
else
  m = 0;
  sz = 0;
  for i = 1 : size( aux, 1 )
    tmp = length(aux{i,1})-1;
    m = m + tmp;
    sz = sz + tmp^2;
  end
  M = spalloc( m, m, sz );
  pos = 0;
  W_pos = 0;
  for i = 1 : size( aux, 1 )
    tmp = length(aux{i,1})-1;
    W_tmp = size(aux{i,3},1);
    I = pos+(1:tmp);
    I_W = W_pos + (1:W_tmp);
    M( I, I ) = lowrankschur6( W(I_W,I_W), aux( i, : ) );
	W_pos = W_pos + W_tmp;
    pos = pos + tmp;
  end
end
