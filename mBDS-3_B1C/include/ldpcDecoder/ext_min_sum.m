function Ls = ext_min_sum(L1, L2)
NM_EMS = 4;
N_GF = 6;           % Number of GF(q) bits
Q_GF = 2^N_GF; 
    if isempty(L1)
        Ls = L2;
        return;
    end
    [~, idx1] = sort(L1);
    [~, idx2] = sort(L2);
    maxL = L1(idx1(NM_EMS)) + L2(idx2(NM_EMS));
    Ls = maxL * ones(1, Q_GF, 'single');

    for i = 1:NM_EMS
        for j = 1:NM_EMS
            xor_idx = bitxor(idx1(i)-1, idx2(j)-1) + 1;  % Convert through zero-based GF indexing.
            if L1(idx1(i)) + L2(idx2(j)) < Ls(xor_idx)
               Ls(xor_idx) = L1(idx1(i)) + L2(idx2(j));
            end
        end
    end
end
