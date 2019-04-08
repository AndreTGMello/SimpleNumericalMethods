function solve_gauss!(A::AbstractArray, b::AbstractArray; verbose::Bool=false)::AbstractArray
    m,n = size(A)
    P = zeros(Int64, m)
    if verbose
        println("----- Iteration 0 -----\n")
        println("--- A ---")
        pprintlnmat(A)
        println("--- b ---")
        pprintlnmat(b)
    end
    for k = 1:m-1
        """
            Exchange rows in case the value at the diagonal is equal to zero.
            Do this for rows in A and b.
        """
        if A[k,k] == 0
            q = k+1
            while A[q,k] == 0 && i < m
                q += 1
            end
            """
                If no other element different from zero can be found, the matrix A
                    is LD and the system can't be solved.
            """
            if A[q,k] == 0
                return [0]
            end
            for i = 1:m
                aux = A[k,i]
                A[k,i] = A[q,i]
                A[q,i] = aux
            end
            aux = b[k]
            b[k] = b[q]
            b[q] = aux
        end
        """
            Loop used for setting all elements in the lines directly below the kth to zero.
            First, use this position (A[i,k]) to store the multiplier A[i,k]/A[k,k],
                then use it for subtracting the kth line from the ith.
            Finaly set A[i,k] to zero.
        """
        for i = k+1:m
            A[i,k] /= A[k,k]
            for j = k+1:m
                A[i,j] -= A[k,j]*A[i,k]
            end
            b[i] -= b[k]*A[i,k]
            A[i,k] = 0
        end
        if verbose
            println("----- Iteration $k -----\n")
            println("--- A ---")
            pprintlnmat(A)
            println("--- b ---")
            pprintlnmat(b)
        end
    end
    """ Solve the equation using A, now triangular. """
    sol = zeros(Float64, (m,1))
    for i = m:-1:1
        sol[i] = b[i]
        for j = m:-1:i+1
            sol[i] -= A[i,j]*sol[j]
        end
        sol[i] /= A[i,i]
    end
    if verbose
        println("----- Finish -----\n")
        println("--- A ---")
        pprintlnmat(A)
        println("--- b ---")
        pprintlnmat(b)
        println("--- x ---")
        pprintlnmat(sol)
    end
    return sol
end

""" --- Tests for Gauss Elimination -Begind --- """
A = Float64[[1 -1 3];[-2 2 4];[1 2 1]]
b = Float64[1 -2 1]'
A\b
res_gauss = solve_gauss!(A,b,verbose=true)
A
b
res = A\b
isapprox(res, res_gauss)

A = rand(-10.0:30.0, (10,10))
b = rand(-1.0:10.0, 10)
A\b
res_gauss = solve_gauss!(A,b)
A
b
res = A\b
isapprox(res, res_gauss)

""" --- Tests for Gauss Elimination -End --- """
