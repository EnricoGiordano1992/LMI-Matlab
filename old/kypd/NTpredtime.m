function t = NTpredtime( stats )

tb = stats.FunctionTable;

for i = 1 : length( tb )
  if strcmp( tb( i ).FunctionName, 'NTpred' )
    t = tb( i ).TotalTime;
    return
  end
end
