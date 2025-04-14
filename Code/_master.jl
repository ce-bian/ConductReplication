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

# Figure 9 
include("CES_extended_3rd_mkt_ad_valorem/main/C2N2_example1.jl") 

# Figure 10
include("CES_extended_3rd_mkt_ad_valorem/main/C4N4_example1.jl")


# main DGP and data output for C2N1 and C3N1
include("CES_MC/main/new_DGP_output.jl")

 
# estimation for varying M with C=2, N=1
include("CES_extended_3rd_mkt_ad_valorem/w_simDT/main/est_simDT_varyM_C2N1.jl")

# optimal trade policy and welfare analysis for C=2, N=1
include("CES_extended_3rd_mkt_ad_valorem/w_simDT/main/policy_table_C2N1.jl")

# Table 1
include("CES_extended_3rd_mkt_ad_valorem/w_simDT/main/parameters_texout_C2N1.jl")

# Figure 5-7
include("CES_extended_3rd_mkt_ad_valorem/w_simDT/main/output_C2N1.jl")

# estimation for varying M with C=3, N=1
include("CES_extended_3rd_mkt_ad_valorem/w_simDT/main/est_simDT_varyM_C3N1.jl")

# optimal trade policy and welfare analysis for C=3, N=1
include("CES_extended_3rd_mkt_ad_valorem/w_simDT/main/policy_table_C3N1.jl")


# Figure 8 
include("CES_extended_3rd_mkt_ad_valorem/w_simDT/main/output_C3N1.jl")