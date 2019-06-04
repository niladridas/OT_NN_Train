function OT_update(Initsamples, Finalsamples, weight_gen)
    # Initsamples: state dim x no. of samples
    M = size(Initsamples)[2]
    w = weight_gen/sum(weight_gen) # weights normalized
    # Calculate the M x M distance matrix
    D = zeros(M,M)
    for i = 2:M # Iterating over row
        for j = (i+1):M # Iterating over column
            D[i,j] = norm(Initsamples[:,i]-Initsamples[:,j]))
            D[j,i] = D[i,j]
        end
    end
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
    Finalsamples = zeros(M)
    for i4= 1:M
        Finalsamples[:,i4] = sum(p_mat[i4,j4]*Initsamples[:,j4] for j4=1:M)
    end
end
