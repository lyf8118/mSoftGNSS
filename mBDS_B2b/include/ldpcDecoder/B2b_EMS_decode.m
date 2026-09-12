function code = B2b_EMS_decode(L)  
    init_table();
    [H_idx,H_ele,H]=getH_information();
    N_GF = 6;           % Number of GF(q) bits
    Q_GF = 2^N_GF;      % Number of GF(q) elements
    MAX_ITER = 15;      % Max number of iterations
    % NM_EMS = 4;         % LLR truncation size of EMS

    n= 162;
    [~, min_indices] = min(L, [], 2); 
    code = min_indices - 1;          
    
    code = code';   

    % Tanner graph edges
    [ie, je, he, ne] = graph_edge(H_idx, H_ele);
    V2C = zeros(ne, Q_GF, 'single');
    C2V = zeros(ne, Q_GF, 'single');


    for i = 1:ne
        V2C(i,:) = permute_V2C(he(i), L(je(i),:));
    end

    for iter = 1:MAX_ITER  

        c = gf(code, 6);        
        Ho = gf(H.', 6); 

        s = c * Ho;

        if all(s == 0)
           
             break;
        end
        % Update check nodes
        for i = 1:ne
            Ls = [];
            for j = 1:ne
                if ie(i) == ie(j) && i ~= j
                    Ls = ext_min_sum(Ls, V2C(j,:));
                end
            end
            Ls = Ls - min(Ls);
            C2V(i,:) = permute_C2V(he(i), Ls);
        end
         
   
        % Update variable nodes
        for i = 1:ne
            Ls = L(je(i),:);
            for j = 1:ne
                if je(i) == je(j) && i ~= j
                    Ls = Ls + C2V(j,:);
                end
            end
            Ls = Ls - min(Ls);
            V2C(i,:) = permute_V2C(he(i), Ls);
        end

        % Update LLR and GF(q) codes
        for i = 1:n
            for j = 1:ne
                if i == je(j)
                    L(i,:) = L(i,:) + C2V(j,:);
                end
            end
            L(i,:) = L(i,:) - min(L(i,:));
            [~, idx] = min(L(i,:));
            code(i) = idx - 1;
        end

    end

