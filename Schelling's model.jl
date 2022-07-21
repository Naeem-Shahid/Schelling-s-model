using ImageFiltering
function segregate(ϕ, p)
    ### Get info for each site
    x = ϕ .== 1; y = ϕ .== 2; vac = ϕ .== 0        # Boolian tag x, y and vacant
    ker = [1.0 1.0 1.0; 1.0 0.0 1.0; 1.0 1.0 1.0]  # Kernel
    xnb = imfilter(x, ker, "circular")             # Compute neighbors
    ynb = imfilter(y, ker, "circular")
    fx = xnb ./ (xnb + ynb)                        # Fraction of neighbors for x and y 
    fy = ynb ./ (xnb + ynb)
    
    ### At this point the new lattice will be constructed with each site  
    ψ = ifelse.(x, fx, fy)                         # Matrix with neighbors of the same kind
    ψ[vac] .= NaN
    unhappy_lst = findall(ψ .< p)
    vac_lst = findall(vac)
    n_vac = sum(vac)
    for site in unhappy_lst                        # Fill sites accordingly
        n = rand(1:n_vac)
        n_site = vac_lst[n]
        ϕ[n_site] = ϕ[site]
        ϕ[site] = 0
        vac_lst[n] = site
    end
end