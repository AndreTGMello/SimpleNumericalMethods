function cholesky_decomposition!(A::AbstractArray)
    m,n = size(A)
    for k = 1:n
        for j = 1:k-1
            A[k,k] -= A[k,j]^2
        end
        A[k,k] = A[k,k]^(1/2)
        for j = k+1:n
            for i = 1:k-1
                A[j,k] -= A[j,i]*A[k,i]
            end
            A[j,k] /= A[k,k]
        end
        for j = k+1:n
            A[k,j] = 0.0
        end
    end
end

""" --- Tests for Cholesky Decomposition -Begin --- """

A = Float64[[4 1 1]; [1 2 -1]; [1 -1 3]]
A_ = copy(A)
lal.det(A)
lal.isposdef(A)
cholesky_decomposition!(A)
A
isapprox(A*A',A_)

""" --- Tests for Cholesky Decomposition -End --- """
