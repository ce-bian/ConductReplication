# ConductReplication
**work in progress, do not use the code yet!**

Replication package for "Inferring Conduct to Guide Strategic Trade Policy" by Eileen (Ce) Bian, Keith Head and Scott Orr

This project is entirely coded in Julia. **Project.toml** contains all required packages and dependencies based on Julia v1.10.0 (the code may or may not work with other versions).

You can install Julia v1.10.0 from Terminal using the following command:
```bash
juliaup add 1.10.0
```
You can also make 1.10.0 the default version:
```bash
juliaup default 1.10.0
```

To install the dependencies, you can either (1) open the folder in VSCode and run the following command in the Julia REPL (this is also included in **install_packages.jl**):
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate() 
```
or (2) open a terminal and navigate to the folder and run the following command:
```bash
julia --project=/Users/cbian/Dropbox/MyProject -e 'using Pkg; Pkg.instantiate()'
```

**_master.jl** is the main file that runs the entire project. It will replicate all tables and figures in the paper. It will take a few days to run (depending on your computer), so you can comment out the parts you don't need. You can also run each file individually by running them in the Julia REPL or in a terminal.

Note that two files are not included in this repository (but will be in the journal replication package) due to size constraints (> 100MB), you can download them from the following links: 
- Data/eqbaB_eqbaC_GP_C2_N1_M1000.jld2: the simulated equilibrium data for C=2 and N=1 
- Data/eqbaB_eqbaC_GP_C3_N1_M1000.jld2: the simulated equilibrium data for C=3 and N=1 

Contact [Eileen Bian](mailto:ce.bian@sauder.ubc.ca) if you have any questions or issues with the code.