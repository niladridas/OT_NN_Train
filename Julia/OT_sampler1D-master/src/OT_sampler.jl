include("../Utils/hist2D.jl")
function OT_sampler(plot_flag,savedata_flag)
    #### Test Examples
    ## Example 1: Beta distribution
    alp = 2
    bet = 5
    B = gamma(alp)*gamma(bet)/gamma(alp+bet)
    f1(x) = x^(alp-1)*(1-x)^(bet-1)/B
    range1 = [0 1]
    M1 = 100
    ## Example 2: Normal Distribution
    mu = 0
    sig2 = 1
    f2(x) = (1/sqrt(2*pi*sig2))*exp((-1/2)*((x-mu)^2)/sig2)
    range2 = [-6 6]
    M2 = 200
    ## Example 3: Cauchy Distribution
    gam = 0.5
    x_0 = 0
    f3(x) = (1/(pi*gam))*(gam^2/((x-x_0)^2+gam^2))
    range3 = [-4 4]
    M3 = 100
    ## Example 4: Bimodal Gaussian distribution
    x_01 = -10
    sigma21 = 2
    x_02 = 10
    sigma22 = 2
    p1 = 0.1
    f4(x) = p1*(1/sqrt(2*pi*sigma21))*exp((-1/2)*((x-x_01)^2)/sigma21) + (1-p1)*(1/sqrt(2*pi*sigma22))*exp((-1/2)*((x-x_02)^2)/sigma22)
    range4 = [-20 20]
    M4 = 200
    ## Example 5: Von Mises Distribution
    k = 8
    I0_k = besseli(0,k)
    f5(x) = exp(k*cos(x))/(2*pi*I0_k)
    range5 = [-1*pi pi]
    M5 = 100
    ##-------------------------------------------------------------------------
    f = Function[f1,f2,f3,f4,f5]
    Mlist = [M1,M2,M3,M4,M5]
    range = [range1,range2,range3,range4,range5]
    for i = 1:length(range)
        # First generate uniform samples from [a,b]
        rangei = range[i]
        a = rangei[1]
        b = rangei[2]
        M = Mlist[i]
        uni_samples = [a+(b-a)*i0/(M+1) for i0 in 1:M] # inline for function
        fi = f[i]
        weight_gen = map(fi,uni_samples)
        w = weight_gen/sum(weight_gen) # weights normalized
        D = abs.(repeat(uni_samples,1,M)-repeat(uni_samples',M,1))
        D_vec = reshape(D,M*M,1) # read along the column
        A1 = zeros(M,M*M)
        A2 = zeros(M,M*M)
        for i1 = 1:M
            A1[i1,M*(i1-1)+1:M*i1] = ones(1,M)
            for j1 = 1:M
                A2[1+j1-1,M*(i1-1)+1+j1-1] = 1
            end
        end
        b2 = ones(M)
        myModel = Model(with_optimizer(GLPK.Optimizer)) #MODEL CONSTRUCTION
        @variable(myModel, p[i2=1:(M*M)] >= 0) # Models x >=0
        @objective(myModel, Min, sum(D_vec[j2]*p[j2] for j2=1:M*M)) # Sets the objective to be minimized. For maximization use Max
        for i3=1:M # for all rows do the following
            @constraint(myModel, sum(A1[i3,j3]*p[j3] for j3=1:M*M) == w[i3]*M)
            @constraint(myModel, sum(A2[i3,j3]*p[j3] for j3=1:M*M) == b2[i3])
        end
        status = @time optimize!(myModel) # solves the model
        p_mat = reshape(value.(p),M,M)
        obj_val = getobjectivevalue(myModel)
        ##
        x_a = zeros(M)
        for i4= 1:M
            x_a[i4] = sum(p_mat[i4,j4]*uni_samples[j4] for j4=1:M)
        end
        if plot_flag == 1
            plotly()
            Plots.histogram(x_a, normed=true,
                legend=:none,
                titlefont=font(9),
                color=:red,
                bins=60)
            Plots.plot!(fi,a,b,color = "blue",show = true)#,reuse = false,show = true) # PGF Plotting
        end
        if savedata_flag == 1
            # Draw a 2D histogram
            (bins , xpoints) = hist2D(x_a,80)
            xmidpts = (xpoints[1:(end-1),1]+xpoints[2:end,1])./2
            data2d  = [xmidpts'; bins']'
            binlength = (xpoints[2]-xpoints[1])
            binnorm = binlength*sum(data2d[:,2])
            data2d[:,2] = data2d[:,2]./binnorm
            str1 = @sprintf("../data/out%i.dat",Int(i))
            writedlm(str1, data2d, " ")
        end
    end
    return 1
end
