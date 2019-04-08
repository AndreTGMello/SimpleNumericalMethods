using Printf
function pprintlnmat(A::AbstractArray; tabs::String="\t\t")
    m,n = size(A)
    print("\n")
    for i = 1:m
        print("|",tabs)
        for j = 1:n
            print("$(@sprintf("%.3f", A[i,j]))",tabs)
        end
        print("|\n")
    end
    print("\n")
end

a = [-1 2 1 3]
pprintlnmat(a)
