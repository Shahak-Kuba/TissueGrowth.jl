using QuadGK
using Roots


function equidistant_points_on_polar_curve(x_function, y_function, num_points)

    function numerical_derivative(f, θ, h=1e-7)
        return (f(θ + h) - f(θ - h)) / (2h)
    end

    # Define the integrand for the arc length in polar coordinates
    integrand = θ -> sqrt(numerical_derivative(x_function, θ)^2 + numerical_derivative(y_function, θ)^2)

    function arc_length(θ)
        result, _ = quadgk(integrand, 0, θ)
        return result
    end

    # Equally spaced points along the polar curve in terms of arc length
    L, _ = quadgk(integrand, 0, 2π)  # Total length of the curve
    Δl = L / (num_points)
    Δθ = 2π / (num_points)

    theta_points = Float64[0.0]
    current_length = 0.0

   rootsFunc = (θ,curr_length) -> arc_length(θ) - (curr_length + Δl)


    for i in 1:num_points - 1
        θ = find_zero(θ->rootsFunc(θ,current_length), (theta_points[i], theta_points[i] + 2*Δθ))
        push!(theta_points, θ)
        current_length = arc_length(θ)
    end

    # Get equally spaced θ values along the polar curve
    θ_values = theta_points

    # Calculate corresponding (x, y) values
    x_values = x_function.(θ_values)
    y_values = y_function.(θ_values)

    return hcat(x_values, y_values)
end

"""
    lineIntersection(rₘ₁, rₗ, rₘ₂, rᵣ)

Implementation of geometric solution for moving a cell in the normal direction. Calculate the intersection points of two lines defined by points `rₘ₁`, `rₗ` and `rₘ₂`, `rᵣ`.

This function computes the intersection point of each pair of lines. If the lines are parallel (the determinant is zero), it returns the midpoint of `rₘ₁` and `rₘ₂`. Otherwise, it calculates the intersection point using the parameters `u` and `t`. If the intersection lies within the segments, the intersection point is returned; otherwise, the midpoint is returned.

# Arguments
- `rₘ₁`: Starting point of the first line segment.
- `rₗ`: Ending point of the first line segment.
- `rₘ₂`: Starting point of the second line segment.
- `rᵣ`: Ending point of the second line segment.

# Returns
A vector of intersection points for each pair of line segments.
"""
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

"""
    Ω(p)

Calculate the area of a polygon defined by points in `p`.

This function computes the area using the shoelace formula. The polygon is defined by a set of points `p`, and the function iterates through these points to calculate the area.

# Arguments
- `p`: A matrix where each row represents a point of the polygon in 2D space.

# Returns
The absolute area of the polygon.
"""
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

"""
    ωκ(rᵢ₋₁, rᵢ, rᵢ₊₁)

Calculate the areas of triangles formed by consecutive triplets of points.

This function takes three points `rᵢ₋₁`, `rᵢ`, and `rᵢ₊₁` and calculates the area of the triangle formed by these points. It is typically used in a loop over a sequence of points to calculate the area of each triangle formed by consecutive triplets.

# Arguments
- `rᵢ₋₁`: The first point of the triangle.
- `rᵢ`: The second point of the triangle.
- `rᵢ₊₁`: The third point of the triangle.

# Returns
A vector of areas for each triangle formed by the input points.
"""
function ωκ(rᵢ₋₁, rᵢ, rᵢ₊₁)
    triVector = [rᵢ₋₁ rᵢ rᵢ₊₁]
    A = zeros(size(triVector,1))
    for ii in axes(triVector,1)
        A[ii] = Ω(reshape(triVector[ii,:],(2,3))')
    end
    return A
end
