;;** AIM: return 16th, 50th and 84th percentile of a distribution

FUNCTION VW_PERCENTILE, X
N = n_elements(X)

X2 = X[sort(X)]
PER16 = X2[N*0.16]
PER50 = X2[N*0.50]
PER84 = X2[N*0.84]

result = [PER16,PER50,PER84]

RETURN, result
END
