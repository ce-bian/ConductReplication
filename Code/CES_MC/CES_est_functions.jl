
function est_demand_OLS(DT::DataFrame)
    """
    Simple OLS regression - Demand estimation
    Input:
        - DT: a dataframe
    Output:
        - estimated σ    
    """
    # OLS regression
    X = [ones(size(DT,1)) log.(DT[!, :P])]
    Y = log.(DT[!, :R_normalized])
    coef = (X'*X)\(X'*Y)
    return 1-coef[2]
end



function est_demand_IV(DT::DataFrame)
    """
    IV regression, with τ as instruments - Demand estimation
    Input:
        - DT: a dataframe
    Output:
        - estimated σ
    """
    # 2SLS regression with τ as instruments
    X = [log.(DT[!, :P]) ones(size(DT,1))]
    Y = log.(DT[!, :R_normalized])
    Z = [ones(size(DT,1)) log.(DT[!, :τ])]
    Q = Z'*X
    W = inv(Z'*Z)
    coef = inv(Q'*W*Q)*Q'*W*(Z'*Y)
    return 1-coef[1]
end



function est_supply_OLS(DT::DataFrame)
    """
    Simple OLS regression - Supply estimation
    Input:
        - DT: a dataframe
    Output:
        - estimated σ    
    """
    # OLS regression
    # X is a matrix with 1 and log(τ)
    X = [ones(size(DT,1)) log.(DT[!, :τ])]
    Y = log.(DT[!, :R_normalized])
    coef = (X'*X)\(X'*Y)
    return 1-coef[2]
end



function est_supply_GMM_twostep(GP::global_param, μR_icm_R::F1, μR_0mm_R::F2,
    DT::DataFrame; market_list_input=missing, seTRUE=true) where {F1 <: Function, F2 <: Function}
    """
    Nonlinear GMM two-step - Supply estimation
    moment condition: E[ln(τ) v_true] = 0
    Input:
       - GP: global parameters
       - DT: a dataframe
       - W: weighted matrix (default is identity matrix)
       - μR_icm_R: based on R, μR_icm_R is either μB_icm_R or μC_icm_R
       - μR_0mm_R: based on R, μR_0mm_R is either μB_0mm_R or μC_0mm_R
       - market_list_input: list of markets to be used for estimation     
       - seTRUE: whether to calculate standard error  
    Output:
        - σ_optimal: estimated σ   
        - se: standard error 
    """
    σ_true = GP.σ
    C = GP.C
    N = GP.N
    M = GP.M
    if ismissing(market_list_input)
        market_list = 1:C+M
    else
        market_list = market_list_input
    end

    n = size(DT, 1) # Number of observations
    Y = log.(DT[!, :R_normalized])
    τ = DT[!, :τ]
    τ1 = hcat(ones(n), log.(τ))  # Use hcat instead of creating arrays in a loop

    RS_long = DT[!, :RS]
    RS = reshape(RS_long, C, N, length(market_list))  # 3-D array; Reshape directly instead of loop assignment
    RS0_long = DT[!, :RS0]

    # set up moment condition
    g = (y, τ1, τ, μ, μ0, σ) -> τ1 * (y - σ[2] - (1 - σ[1]) * log(τ) - (1 - σ[1]) * (log(μ)-log(μ0)))

    function sample_moment(RS, RS0_long, σ_est)
        μ_array = [μR_icm_R(c, i, RS[:, :, m], σ_est[1]) for c in 1:C, i in 1:N, m in 1:length(market_list)] 
        μ0_array = [μR_0mm_R(RS0_long[i], σ_est[1]) for i in 1:n] 
        # μ = μ_array
        V = mean(g(Y[i], τ1[i, :], τ[i], μ_array[i], μ0_array[i], σ_est)' for i in 1:n)  # Use mean instead of sum/n
        return V
    end

    function sample_moment_matrix(RS, RS0_long,σ_est)
        μ_array = [μR_icm_R(c, i, RS[:, :, m], σ_est[1]) for c in 1:C, i in 1:N, m in 1:length(market_list)]
        μ0_array = [μR_0mm_R(RS0_long[i], σ_est[1]) for i in 1:n]
        #μ = μ_array
        V = hcat((g(Y[i], τ1[i, :], τ[i], μ_array[i],μ0_array[i], σ_est) for i in 1:n)...)  # Use generator instead of temporary array
        return V
    end


    function obj(ι_est,RS, RS0_long, W)
        """
        objective function for nonlinear-GMM
        """
        σ_est = ι_est
        σ_est[1] = 1+exp(ι_est[1]) # ensure σ>1 -> μ>0
        V = sample_moment(RS, RS0_long, σ_est)
        return V * W * V'
    end
    # Initial guess for σ_est = σ_true
    ι_est_init = [log(σ_true-1) 0.0]
    # Initial weighting matrix is the identity matrix
    W1 = I(2)
    options = Optim.Options(g_tol=1e-12, f_tol=1e-12, x_tol=1e-12)
    result = optimize(ι -> obj(ι, RS, RS0_long, W1), ι_est_init, LBFGS(), options, autodiff=:forward)
    ι_est = Optim.minimizer(result)
    σ_est = ι_est
    σ_est[1] = 1+exp(ι_est[1]) 
    # Step 2: Compute Optimal Weighting Matrix and Re-estimate
    # Create a matrix of the sample moments at the initial estimate
    V = sample_moment_matrix(RS, RS0_long,σ_est)
    # Compute the optimal weighting matrix
    W2 = inv(V * V' / n)
    # Re-estimate the parameters using the optimal weighting matrix
    result = optimize(ι -> obj(ι,RS, RS0_long, W2), ι_est_init, LBFGS(),options, autodiff=:forward) # investigate later, why here ι_est doesn't work (sometimes it gives σ_optimal = 1.0 which means ι is negative infinity)
    ι_optimal = Optim.minimizer(result)
    σ_optimal = ι_optimal
    σ_optimal[1] = 1+exp(ι_optimal[1])

    if seTRUE
        # compute standard error
        # Get asyvar, we compute the gradient of g with respect to theta
        # Jacobian will yield dg(W,θ)/dθ'
        μ_array = [μR_icm_R(c, i, RS[:, :, m], σ_optimal[1]) for c in 1:C, i in 1:N, m in 1:length(market_list)]
        μ0_array = [μR_0mm_R(RS0_long[i], σ_optimal[1]) for i in 1:n]
        # μ = μ_array
        v = map(i -> ForwardDiff.jacobian(θ -> g(Y[i], τ1[i, :], τ[i], μ_array[i], μ0_array[i], θ), σ_optimal), 1:n)
        dg = sum([v[i] for i in 1:n])

        V = sample_moment_matrix(RS, RS0_long, σ_optimal)
        W = inv(V * V')
        
        avar = inv(dg' * W * dg)
        se = sqrt.(diag(avar))
        return σ_optimal, se
    else
        return σ_optimal, [missing, missing]
    end
end



function est_supply_demand_GMM_twostep(GP::global_param, μR_icm_R::F1, μR_0mm_R::F2, DT::DataFrame; market_list_input=missing, seTRUE= true) where {F1 <: Function, F2 <: Function}
    """
    Two-step Nonlinear GMM - Supply and demand, using LBFGS()
    two moment conditions: (1) E[ln(τ) v_true] = 0 (2) E[ln(τ) ξ_true] = 0
    Input:
        - GP: global parameters
        - μR_icm_R: μB_icm_R (Bertrand) or μC_icm_R (Cournot)
        - μR_0mm_R: μB_0mm_R (Bertrand) or μC_0mm_R (Cournot)
        - DT: a dataframe
        - market_list_input: list of markets to be used for estimation
        - seTRUE: whether to calculate standard error
    Output:
        - σ_optimal: estimated σ
        - se: standard error
    """
    σ_true = GP.σ
    C = GP.C
    N = GP.N
    M = GP.M

    market_list = market_list_input
    # check whether we have specific market_list
    if ismissing(market_list_input)
        market_list = 1:C+M
    else
        market_list = market_list_input
    end

    n = size(DT, 1)
    Y = log.(DT[!, :R_normalized])
    P = DT[!, :P]
    τ = DT[!, :τ]
    τ1 = hcat(ones(n), log.(τ))  # Use hcat instead of creating arrays in a loop

    RS_long = DT[!, :RS]
    RS = reshape(RS_long, C, N, length(market_list))  # Reshape directly instead of loop assignment
    RS0_long = DT[!, :RS0]

    function g(y, p, τ1, τ, μ, μ0, σ)
        """
        # set up moment condition
        """
        return [τ1 * (y - σ[2] - (1 - σ[1]) * log(p));
                τ1 * (y - σ[3] - (1 - σ[1]) * log(τ) - (1 - σ[1]) * (log(μ)-log(μ0)))]'
    end

    function sample_moment(RS, RS0_long,σ_est)
        """
        calculate the sample moment based on σ_est, return a vector
        """
        μ_array = [μR_icm_R(c, i, RS[:, :, m], σ_est[1]) for c in 1:C, i in 1:N, m in 1:length(market_list)]
        μ0_array = [μR_0mm_R(RS0_long[i], σ_est[1]) for i in 1:n]
        # μ = μ_array
        V = mean(g(Y[i], P[i], τ1[i, :], τ[i], μ_array[i], μ0_array[i], σ_est) for i in 1:n)  # Use mean instead of sum/n
        return V
    end
    function sample_moment_matrix(RS, RS0_long,σ_est)
        """
        calculate the sample moment based on σ_est, return a matrix
        """
        μ_array = [μR_icm_R(c, i, RS[:, :, m], σ_est[1]) for c in 1:C, i in 1:N, m in 1:length(market_list)]
        μ0_array = [μR_0mm_R(RS0_long[i], σ_est[1]) for i in 1:n]
        # μ = μ_array
        V = hcat((g(Y[i], P[i], τ1[i, :], τ[i], μ_array[i], μ0_array[i], σ_est)' for i in 1:n)...)  # Use generator instead of temporary array
        return V
    end

    function obj(ι_est, RS, RS0_long, W)
        """
        objective function for nonlinear-GMM
        """
        σ_est = ι_est
        σ_est[1] = 1+exp(ι_est[1]) # ensure σ>1 -> μ>0
        V = sample_moment(RS, RS0_long, σ_est)
        return V * W * V'
    end

    ι_est_init = [log(σ_true-1) 0.0 0.0]
    W1 = I(4)
    options = Optim.Options(g_tol=1e-12, f_tol=1e-12, x_tol=1e-12)
    result = optimize(ι -> obj(ι, RS, RS0_long, W1), ι_est_init, LBFGS(), options, autodiff=:forward)
    ι_est = Optim.minimizer(result)
    σ_est = ι_est
    σ_est[1] = 1+exp(ι_est[1]) 
    
    # second step
    V = sample_moment_matrix(RS, RS0_long, σ_est)
    W2 = inv(V * V' / n)
    result = optimize(ι -> obj(ι, RS, RS0_long, W2), ι_est_init, LBFGS(),options, autodiff=:forward) # investigate later, why here ι_est doesn't work (sometimes it does, sometimes it doesn't, give σ_optimal = 1.0 which means ι is negative infinity)
    ι_optimal = Optim.minimizer(result)
    σ_optimal = ι_optimal
    σ_optimal[1] = 1+exp(ι_optimal[1])
    

    if seTRUE
        # compute standard error
        # Get asyvar, we compute the gradient of g with respect to theta
        # Jacobian will yield dg(W,θ)/dθ'
        μ_array = [μR_icm_R(c, i, RS[:, :, m], σ_optimal[1]) for c in 1:C, i in 1:N, m in 1:length(market_list)]
        μ0_array = [μR_0mm_R(RS0_long[i], σ_optimal[1]) for i in 1:n]
        # μ = μ_array
        v = map(i -> ForwardDiff.jacobian(θ -> g(Y[i], P[i], τ1[i, :], τ[i], μ_array[i], μ0_array[i], θ), σ_optimal), 1:n)
        dg = sum([v[i] for i in 1:n])

        V = sample_moment_matrix(RS, RS0_long, σ_optimal)
        W = inv(V * V')
        
        avar = inv(dg' * W * dg)
        se = sqrt.(diag(avar))
        return σ_optimal, se
    else
        return σ_optimal, [missing, missing, missing]
    end
end
