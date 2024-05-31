function printInfo(simNum,simTotal,btype,N,kₛ,η,kf,M,D)
    println(@sprintf "----------------------------- Simulation %d/%d Complete -----------------------------" simNum simTotal)
    println(@sprintf "Boundary Type: %s, Cell count: %d, Springs per cell: %d" btype N Int(M/N))
    println(@sprintf "kₛ¹: %.5f, η¹: %.5f, kf¹: %.5f, Diffusivity: %.5f" kₛ η kf D)
    println(@sprintf "-----------------------------------------------------------------------------------")
end

# Helper function to convert ElasticMatrix and ElasticVector to regular Matrix and Vector
function convert_to_regular_matrix_vector(sim_results::SimResults_t)
    u_x_positions = [el[:, 1].data for el in sim_results.u]
    u_y_positions = [el[:, 2].data for el in sim_results.u]
    ∑F_vectors = [el.data for el in sim_results.∑F]
    Density_matrices = [el.data for el in sim_results.Density]
    ψ_vectors = [el.data for el in sim_results.ψ]
    Κ_matrices = [el.data for el in sim_results.Κ]
    
    return (u_x_positions, u_y_positions, ∑F_vectors, Density_matrices, ψ_vectors, Κ_matrices)
end

# Function to convert a vector of vectors to a DataFrame column of vectors
function vector_of_vectors_to_column(vov::Vector{Vector{T}}) where T
    return vov
end

# Function to export the data to CSV
function export_data_to_csv(sim_results::SimResults_t, filename::String)
    u_x_positions, u_y_positions, ∑F_vectors, Density_matrices, ψ_vectors, Κ_matrices = convert_to_regular_matrix_vector(sim_results)

    # Convert vectors of vectors to DataFrame columns
    df = DataFrame(
        t = sim_results.t,
        pos_x = vector_of_vectors_to_column(u_x_positions),
        pos_y = vector_of_vectors_to_column(u_y_positions),
        Density = vector_of_vectors_to_column(Density_matrices),
        Stress = vector_of_vectors_to_column(ψ_vectors),
        Force_Sum = vector_of_vectors_to_column(∑F_vectors),
        Velocity = vector_of_vectors_to_column(sim_results.Vₙ),
        Void_Area = sim_results.Ω,
        Curvature = vector_of_vectors_to_column(Κ_matrices),
        Cell_Count = sim_results.CellCount
    )
    
    CSV.write(filename, df)
end
