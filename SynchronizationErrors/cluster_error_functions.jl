#! /opt/bin/julia 

using MAT
using LinearAlgebra
using Statistics


# Compute the error corresponding to a cluster
function ClusterError(soly);
    
    meany  = mean(soly,dims=1);
    resta = (soly .- meany).^2;
    raiz = sqrt.(mean(resta,dims=1));
    Err  = mean(raiz,dims=2)[:,1];
    
    return Err[1];
end



# Iterate over each cluster in cluster_indices, computing their errors
function clusynch(coupling, cluster_indices, sim);
    
    Error_dict = Dict()
    
    # Error Synch Clusters
    for (cluster, indices) in cluster_indices
        listsol = [sim[x,:,coupling] for x in indices];
        soly = mapreduce(permutedims, vcat, listsol[:][:]);
        Error_dict[cluster] = ClusterError(soly);
    end

    # Error Synch Global
    solyT  = sim[:,:,coupling];
    Error_dict["ErrT"] = ClusterError(solyT);
    
    return Error_dict
end



# Iterate over the range of couplings, computing the cluster errors.
function compute_errors(cluster_indices, sim)
    Error_dict = clusynch(1, cluster_indices, sim)
    Err_matrix =  collect(values(Error_dict))

    for coupling in 2:len

        Error_dict = clusynch(coupling, cluster_indices, sim)
        Err_matrix = hcat(Err_matrix, collect(values(Error_dict)))

    end
    
    return Float64.(Err_matrix);
end