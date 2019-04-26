
function CdW = CderW(W,Delta,A,input) 
    for i = 1:length(W) 
        Wder = W{i};
        for k = 1:size(Wder,1)
            for j = 1:size(Wder,2)
                if i == 1
                    T = input;
                else
                    T = A{i-1};
                end
                Wder(k,j) = Delta{i}(k,1)*T(j,1);
            end
        end
        CdW{i} = Wder;
    end
end