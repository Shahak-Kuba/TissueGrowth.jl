function lineIntersection(rₘ₁,rₗ,rₘ₂,rᵣ)
    intersect = zeros(size(rₘ₁))
    r = rₗ.-rₘ₁
    s = rᵣ.-rₘ₂
    
    d = r[:,1].*s[:,2] - r[:,2].*s[:,1]
    # performing determinant test in case lines are parallel
    for i in 1:size(d,1)
        if d[i] == 0
            intersect[i,:] = (rₘ₁[i,:] + rₘ₂[i,:])./2
        
        else
            u = ((rₘ₂[i,1] - rₘ₁[i,1])*r[i,2] - (rₘ₂[i,2] - rₘ₁[i,2])*r[i,1])/d[i]
            t = ((rₘ₂[i,1] - rₘ₁[i,1])*s[i,2] - (rₘ₂[i,2] - rₘ₁[i,2])*s[i,1])/d[i]

            if 0≤u≤1 && 0≤t≤1
                #println("Yes these intersect at: ")
                #println(rₘ₁ + t*r)
                intersect[i,:] = (rₘ₁[i,:] + t*r[i,:])
            else
                #println("No these lines dont intersect, midpoint: " )
                #println(((rₘ₁ + rₘ₂)/2))
                intersect[i,:] =  (rₘ₁[i,:] + rₘ₂[i,:])/2
            end
        end
    end
    return intersect
end

function Ω(p)
    A = 0
    for ii in axes(p,1)
        if ii == size(p,1)
            A += (p[ii,1]*p[1,2] -  p[ii,2]*p[1,1])
        else
            A += (p[ii,1]*p[ii+1,2] -  p[ii,2]*p[ii+1,1])
        end
    end
    return abs(A)/2;
end

function ωκ(rᵢ₋₁, rᵢ, rᵢ₊₁)
    triVector = [rᵢ₋₁ rᵢ rᵢ₊₁]
    A = zeros(size(triVector,1))
    for ii in axes(triVector,1)
        A[ii] = Ω(reshape(triVector[ii,:],(2,3))')
    end
    return A
end
