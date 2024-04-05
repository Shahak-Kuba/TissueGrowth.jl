function printInfo(simNum,simTotal,btype,N,kₛ,η,kf,M,D)
    println(@sprintf "----------------------------- Simulation %d/%d Complete -----------------------------" simNum simTotal)
    println(@sprintf "Boundary Type: %s, Cell count: %d, Springs per cell: %d" btype N Int(M/N))
    println(@sprintf "kₛ¹: %.5f, η¹: %.5f, kf¹: %.5f, Diffusivity: %.5f" kₛ η kf D)
    println(@sprintf "-----------------------------------------------------------------------------------")
end

function export_figure(fig, image_name, image_type)
    # Get the current directory
    current_dir = pwd()
    
    # Generate the save directory path
    save_dir = joinpath(current_dir, "figures")

    # Check if save directory exists, create it if it doesn't
    if !isdir(save_dir)
        mkdir(save_dir)
    end
    
    # Generate the timestamp
    timestamp = Dates.now()
    timestamp_str = Dates.format(timestamp, "yyyymmdd")
    
    # Generate the filename
    filename = joinpath(save_dir, "$image_name-$timestamp_str.$image_type")
    
    # Save the figure
    save(filename, fig)
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

