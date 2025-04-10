# Set working directory to ConductReplication
cd(abspath(joinpath(@__DIR__, "..")))

# install the necessary packages
include("install_packages.jl")

# Figure 1
include("CES_3rd_mkt_ad_valorem/main/CES_3rd_mkt_reaction_function_plots.jl")

# Figure 2
include("CES_3rd_mkt_ad_valorem/main/CES_3rd_mkt_example.jl")

# Figure 3
include("CES_extended_3rd_mkt_ad_valorem/main/C2N1_example1.jl")

# Figure 4
include("CES_extended_3rd_mkt_ad_valorem/main/C2N1_example2.jl")



# main DGP and estimation for M = 1000 
include("CES_MC/main/new_DGP_output.jl")


# not run yet =====
# Figure 9 
include("CES_extended_3rd_mkt_ad_valorem/main/C2N2_example1.jl") 
# modification: change Y_m = 1.0 to be consistent with other figures 

# Figure 10 (note that the text in draft is inconsistent)
include("CES_extended_3rd_mkt_ad_valorem/main/C4N4_example1.jl")
