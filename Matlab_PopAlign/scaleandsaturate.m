function df= scaleandsaturate(df,hi,lo)

%  function scaleandsaturate(df,hi,lo)
%  this function scales and saturates a dataset according to the hi and lo
%  thresholds. values < lo get set to 0 and values > hi get set to hi

% df = df ./ std(df')';
df = df ./ max(df')';
df(df < lo) = 0;
df(df > hi) = hi;


