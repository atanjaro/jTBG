"""
    plot_electron_energy( ε1, ε2, ε3, ks1, ks2, ks3, path )    

Generates a plot of the electron band structure along specified cuts of the FBZ.

"""
function plot_electron_energy(el_energies1,el_energies2,el_energies3,ΓM_points,MK_points,KΓ_points,directory_path)
    gr()

    # flatten k-points 
    flat_k1 = permutedims(hcat(ΓM_points...))
    flat_k2 = permutedims(hcat(MK_points...))
    flat_k3 = permutedims(hcat(KΓ_points...))

    # clear plot
    plot(fontfamily="Computer Modern",size=(800, 600),xtickfontsize=14,ytickfontsize=14,xguidefontsize=14,yguidefontsize=14)

    # plot electronic energies from:
    #########################
    ##       Γ to M        ##
    #########################
    plot!(flat_k1[:,1],el_energies1[:,1],linecolor=:red,linewidth=1.5,)
    plot!(flat_k1[:,1],el_energies1[:,2],linecolor=:blue,linewidth=1.5,)
    ########################
    ##       M to K       ##
    ########################
    plot!(flat_k2[:,2].+maximum(flat_k1),el_energies2[:,1],linecolor=:red,linewidth=1.5,)
    plot!(flat_k2[:,2].+maximum(flat_k1),el_energies2[:,2],linecolor=:blue,linewidth=1.5,)
    ########################
    ##       K to Γ       ##
    ########################
    plot!(flat_k3[:,1].+(maximum(flat_k2[:,2].+maximum(flat_k1))),reverse(el_energies3[:,1]),linecolor=:red,linewidth=1.5,)
    plot!(flat_k3[:,1].+(maximum(flat_k2[:,2].+maximum(flat_k1))),reverse(el_energies3[:,2]),linecolor=:blue,linewidth=1.5,)

    vline!([maximum(flat_k1),maximum(flat_k2[:,2].+maximum(flat_k1))], line=:solid, linewidth=0.25, color=:gray)
    plot!(legend=false)
    plot!(xlabel=L"k\:\left[\frac{1}{a}\right]", ylabel=L"ε\:[\textrm{eV}]", xlims=(0,maximum(flat_k3[:,1]).+maximum(flat_k2[:,2].+maximum(flat_k1))), ylims=(-10, 10)) 
    plot!(xticks=([0,maximum(flat_k1),maximum(flat_k2[:,2].+maximum(flat_k1)),maximum(flat_k2[:,2].+maximum(flat_k1).+maximum(flat_k3))],[L"\Gamma",L"M",L"K",L"\Gamma"])) 

    # # Add solid lines for the top and right sides using hline! and vline!
    # hline!([10, 10], line=:solid, linewidth=1.5,color=:black)
    # vline!([maximum(flat_k1), maximum(flat_k2[:,2] .+ maximum(flat_k1)), maximum(flat_k2[:,2] .+ maximum(flat_k1) .+ maximum(flat_k3))], line=:solid, linewidth=1.5,color=:black)

    # Display the final plot
    gui(Plots.plot!())

    # Save the plot to a file
    savefig(directory_path*"graphene_band_structure.pdf") 


    return nothing
end

"""
    plot_phonon_energy( ω1, ω2, ω3, ω4, ω5, ω6, ks1, ks2, ks3, path )    

Generates a plot of the dispersion of all phonon modes along specified cuts of the FBZ.

"""
function plot_phonon_energy(ph_energies1,ph_energies2,ph_energies3,ph_energies4,ph_energies5,ph_energies6,ΓM_points,MK_points,KΓ_points,directory_path)
    gr()

    # flatten k-points 
    flat_k1 = permutedims(hcat(ΓM_points...))
    flat_k2 = permutedims(hcat(MK_points...))
    flat_k3 = permutedims(hcat(KΓ_points...))


    # plot phononic energies from:
    #########################
    ##       Γ to M        ##
    #########################
    # out-of-plane modes
    plot!(flat_k1[:,1],ph_energies1[:,1],linecolor=:red, line=:solid,linewidth=1.5,)
    plot!(flat_k1[:,1],ph_energies1[:,2],linecolor=:red, line=:solid,linewidth=1.5,)
    # in-plane modes
    plot!(flat_k1[:,1],ph_energies4[:,1],linecolor=:black, line=:solid,linewidth=1.5,)
    plot!(flat_k1[:,1],ph_energies4[:,2],linecolor=:black, line=:solid,linewidth=1.5,)
    plot!(flat_k1[:,1],ph_energies4[:,3],linecolor=:black, line=:solid,linewidth=1.5,)
    plot!(flat_k1[:,1],ph_energies4[:,4],linecolor=:black, line=:solid,linewidth=1.5,)
    ########################
    ##       M to K       ##
    ########################
    # out-of-plane modes
    plot!(flat_k2[:,2].+maximum(flat_k1),ph_energies2[:,1],linecolor=:red, line=:solid,linewidth=1.5,)
    plot!(flat_k2[:,2].+maximum(flat_k1),ph_energies2[:,2],linecolor=:red, line=:solid,linewidth=1.5,)
    # in-plane modes
    plot!(flat_k2[:,2].+maximum(flat_k1),ph_energies5[:,1],linecolor=:black, line=:solid,linewidth=1.5,)
    plot!(flat_k2[:,2].+maximum(flat_k1),ph_energies5[:,2],linecolor=:black, line=:solid,linewidth=1.5,)
    plot!(flat_k2[:,2].+maximum(flat_k1),ph_energies5[:,3],linecolor=:black, line=:solid,linewidth=1.5,)
    plot!(flat_k2[:,2].+maximum(flat_k1),ph_energies5[:,4],linecolor=:black, line=:solid,linewidth=1.5,)
    ########################
    ##       K to Γ       ##
    ########################
    # out-of-plane modes
    plot!(flat_k3[:,1].+(maximum(flat_k2[:,2].+maximum(flat_k1))),reverse(ph_energies3[:,1]),linecolor=:red, line=:solid,linewidth=1.5,)
    plot!(flat_k3[:,1].+(maximum(flat_k2[:,2].+maximum(flat_k1))),reverse(ph_energies3[:,2]),linecolor=:red, line=:solid,linewidth=1.5,)
    # in-plane modes
    plot!(flat_k3[:,1].+(maximum(flat_k2[:,2].+maximum(flat_k1))),reverse(ph_energies6[:,1]),linecolor=:black, line=:solid,linewidth=1.5,)
    plot!(flat_k3[:,1].+(maximum(flat_k2[:,2].+maximum(flat_k1))),reverse(ph_energies6[:,2]),linecolor=:black, line=:solid,linewidth=1.5,)
    plot!(flat_k3[:,1].+(maximum(flat_k2[:,2].+maximum(flat_k1))),reverse(ph_energies6[:,3]),linecolor=:black, line=:solid,linewidth=1.5,)
    plot!(flat_k3[:,1].+(maximum(flat_k2[:,2].+maximum(flat_k1))),reverse(ph_energies6[:,4]),linecolor=:black, line=:solid,linewidth=1.5,)

    vline!([maximum(flat_k1),maximum(flat_k2[:,2].+maximum(flat_k1))], line=:solid, linewidth=0.25, color=:gray)
    plot!(legend=false)
    plot!(xlabel=L"k\:\left[\frac{1}{a}\right]", ylabel=L"\Omega\:[\mathrm{eV}]", xlims=(0,maximum(flat_k3[:,1]).+maximum(flat_k2[:,2].+maximum(flat_k1))),ylims=(0,0.2)) 
    plot!(xticks=([0,maximum(flat_k1),maximum(flat_k2[:,2].+maximum(flat_k1)),maximum(flat_k2[:,2].+maximum(flat_k1).+maximum(flat_k3))],[L"\Gamma",L"M",L"K",L"\Gamma"]))

    # out-of-plane mode annotations
    annotate!(1.5, 0.015, L"\mathrm{ZA}", :red, halign=:center, valign=:center)
    annotate!(0.25, 0.1125, L"\mathrm{ZO}", :red, halign=:center, valign=:center)

    # in-plane-mode annotations
    annotate!(0.6, 0.175, L"\mathrm{LO}", :black, halign=:center, valign=:center)
    annotate!(0.6, 0.2, L"\mathrm{TO}", :black, halign=:center, valign=:center)

    annotate!(1.5, 0.1125, L"\mathrm{LA}", :black, halign=:center, valign=:center)
    annotate!(2.6, 0.1125, L"\mathrm{TA}", :black, halign=:center, valign=:center)


    # Display the final plot
    gui(Plots.plot!())


    # Save the plot to a file
    savefig(directory_path*"graphene_phonon_dispersion.pdf") 

    return nothing
end