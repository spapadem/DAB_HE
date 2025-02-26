#function matmupar(A::Array{Float64, 2},x::Array{Float64, 1})
function matmupar(A,x)

	output = SharedArray{Float64}(size(x,1),size(x,2));
    @inbounds @sync @distributed for i in 1:size(A, 1)
        @inbounds for j in 1:size(x, 1)
            output[i] += A[i, j] * x[j]
        end
    end

    return output
end

