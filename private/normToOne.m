function s = normToOne(s)

sSum = sum(s);

n = length(s);
s = s./sSum(ones(n,1));
