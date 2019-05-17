function X = eq_wsamples(X1,W1,W2)
%   eq_wsamples:  Returns equally weighted samples.
% 
%   X = eq_wsamples(X1,W1,W2) 
%       Takes X1 as the prior samples with weights W1.
%       The weights of X1 in the posterior pdf is W2.
%   
%  Dimension of X1 is n (state dimension) x Ns (sample size)
%  Dimension of W1 is 1 x Ns (sample size)
%  Dimension of W2 is 1 x Ns (sample size)
    Ns = size(X1,2);
    n = size(X1,1);
    X = zeros(n,Ns);
    P = opt_map(X1,W1,W2);
    for i = 1:Ns
        for j = 1:Ns
            X(:,i) = X(:,i) + P(i,j)*X1(:,j);
        end
    end
end
% Generates the P matrix from Optimal Map T
function P = opt_map(X1,W1,W2)
    Ns = size(X1,2);
    D = cost_mat(X1);
    D_vec = reshape(D,Ns*Ns,1);
    [A1,A2] = lp_consts(Ns);
    options = optimoptions('linprog','Algorithm','dual-simplex','display', 'off');
    T = linprog(D_vec,[],[],[A2;A1],[W1';W2'],zeros(Ns*Ns,1),[],[],options);
    P = T./repmat(W1',Ns,1); 
    P = reshape(P,[Ns,Ns]);
end
% Cost function
function D = cost_mat(x)
    Ns = size(x,2);
    D = zeros(Ns,Ns);
    for i = 1:Ns
        for j = 1:Ns
            D(i,j) = norm(x(:,i)-x(:,j));
        end
    end
end
% LP constant matrices
function [A1,A2] = lp_consts(Ns)
    A1 = zeros(Ns,Ns*Ns);
    A2 = zeros(Ns,Ns*Ns);
    for i = 1:Ns
        A1(i,Ns*(i-1)+1:Ns*i) = ones(1,Ns);
        for j = 1:Ns
            A2(1+j-1,Ns*(i-1)+1+j-1) = 1;
        end
    end
end