function [OT_samplesX] = OTsink(samples,nSamp,nx,sinkparams,W0,W_post)
    M = cost_mat(samples);
    K=exp(-sinkparams.lambda*M);
    U=K.*M;
    [~,~,u,v]=sinkhornTransport(W0',W_post',K,U,sinkparams.lambda,...
        sinkparams.stoppingCriterion,sinkparams.p_norm,sinkparams.tolerance,sinkparams.maxIter,0);
    T = diag(u(:,1)) * K * diag(v(:,1));
    P = T./repmat(W0,nSamp,1); 
    X = zeros(nx,nSamp);
    for k = 1:nSamp
        for j = 1:nSamp
            X(:,k) = X(:,k) + P(k,j)*samples(:,j);
        end
    end
    OT_samplesX = X;
end
function D = cost_mat(x)
    Ns = size(x,2);
    D = zeros(Ns,Ns);
    for i = 1:Ns
        for j = 1:Ns
            D(i,j) = norm(x(:,i)-x(:,j));
        end
    end
end