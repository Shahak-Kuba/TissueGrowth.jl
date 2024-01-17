using LinearAlgebra
using Printf
using JLD2
import FilePaths

function printInfo(simNum,simTotal,btype,N,kₛ,η,kf,M,D)
    println(@sprintf "----------------------------- Simulation %d/%d Complete -----------------------------" simNum simTotal)
    println(@sprintf "Boundary Type: %s, Cell count: %d, Springs per cell: %d" btype N Int(M/N))
    println(@sprintf "kₛ¹: %.5f, η¹: %.5f, kf¹: %.5f, Diffusivity: %.5f" kₛ η kf D)
    println(@sprintf "-----------------------------------------------------------------------------------")
end

function SaveData(data, SaveName, SaveFolder)
    # Check if the folder exists, create it if it doesn't
    if !isdir(SaveFolder)
        mkdir(SaveFolder)
    end
    
    # Construct the file path
    filepath = joinpath(SaveFolder, "$SaveName.jld2")
    
    # Check if the file already exists
    if isfile(filepath)
        error("File '$SaveName.jld2' already exists in folder '$SaveFolder'. Please choose a different name.")
    else
        # Save the data to the specified file
        save(filepath, "data", data)
    end
end

function LoadData(SaveName, SaveFolder)
    # Construct the file path
    filepath = joinpath(SaveFolder, "$SaveName.jld2")
    
    # Load the data from the specified file
    loaded_data = load(filepath)
    
    # Retrieve the data from the loaded file
    if haskey(loaded_data, "data")
        return loaded_data["data"]
    else
        error("No 'data' key found in the loaded file.")
    end
end

