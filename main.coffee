nm = numeric #;
csvd = $blab.csvd # complex SVD imported in "Ext"

# !!! input matrix must be complex !!!

# Helper functions
size = (A) -> [A.x.length, A.x[0].length]
max = (x) -> Math.max.apply null, x

# Complex pseudoinverse
cpinv = (A) ->
    {U, S, V} = csvd A
    tol = max(size(A))*nm.epsilon*S[0]
    Z = ((if x>tol then 1/x else 0) for x in S)
    V.dot(nm.diag(Z)).dot(U.transjugate())

# Export to other blabs
$blab.cpinv = cpinv

# Test
A = new nm.T([[1, 0, 1],[1, 3, 0]],[[0, 0, 0],[1, 0, 0]]) 

Ap = cpinv(A) #;
# Real part of pseudoinverse
Ap.x
# Imag part of pseudoinverse
Ap.y


