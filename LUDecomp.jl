function LU!(A)
    m,n = size(A)
    P = zeros(Int32, m)
    """ Opperation need to be carried out for each row and column,
    thus the 1:m loop """
    for k = 1:m
        """ First step: update all the rows (fix column, iterate on rows) """
        for j = k:m
            """ Subtract from present value the products from
            all the values already calculated
            that match your row and column """
            for i = 1:k-1
                A[j,k] -= A[i,k]*A[j,i]
            end
        end
        """ Second step: Check if value on the diagonal is different from zero """
        if A[k,k] == 0.0 && k != m
            q = k+1
            while A[q,k] == 0
                q += 1
                if q == m && A[q,k] == 0
                    return "Decomposition failed. Matrix can't be decomposed"
                end
            end
            P[k] = q
            for i = 1:m
                aux = A[q,i]
                A[q,i] = A[k,i]
                A[k,i] = aux
            end
        end
        """ Third step: update all the columns (fix row, iterate on columns);
        Don't update the same value twice!! Begin from k+1.
        Finaly, update all the rows from k+1 up to n """
        for o = k+1:m
            for p = 1:k-1
                A[k,o] -= A[k,p]*A[p,o]
            end
            A[o,k] /= A[k,k]
        end
    end
    return P
end

function solve_PALU(P::AbstractArray, A::AbstractArray, b::AbstractArray)
    m,n = size(A)
    x = zeros(eltype(A), m)
    y = zeros(eltype(A), m)
    for i = 1:n
        if P[i] != 0
            aux = b[P[i]]
            b[P[i]] = b[i]
            b[i] = aux
        end
    end
    """ First solve Ly=b """
    for i = 1:m
        y[i] = b[i]
        for j = 1:i-1
            y[i] -= A[i,j]*y[j]
        end
    end
    """ Then solve Ux=y """
    for i = m:-1:1
        x[i] = y[i]
        println(x[i])
        for j = i+1:m
            x[i] -= A[i,j]*x[j]
        end
        x[i] /= A[i,i]
        println(A[i,i])
        println(x[i])
    end
    return x
end


""" --- Tests for PA=LU Decomposition -Begin --- """
import LinearAlgebra const lal = LinearAlgebra
using BenchmarkTools

A = Float64[[4 3];[6 3]]
@time P = LU!(A)
A
P
b = [1 4]
solve_PALU(P,A,b)
3/2
5/3

A = Float64[[1 -1 3];[-2 2 4];[1 2 1]]
P = LU!(A)
A
P
b = [1 -2 1]
solve_PALU(P,A,b)


A = Float64[[1 2 3];[3 2 1];[5 6 7]]
P = LU!(A)
A
P

A = Float64[[1 -1 3];[1 5 -2];[-2 0 4]]
P = LU!(A)
A
P

A = Float64[[1 2 3];[3 2 11];[5 26 7]]
P = LU!(A)
A
P
linal_fact = lal.lu(A)

A = rand(-10.0:20.0, (3,3))
P = LU!(A)
P
A
linal_fact = lal.lu(A)
linal_fact.factors == A

m = 1000
A = rand(-10.0:50.0, (m,m))
@benchmark P = LU!(A)
@benchmark lal.lu(A)
P
A
linal_fact = lal.lu(A)
linal_fact.factors == A

""" --- Tests for PA=LU Decomposition -End --- """
