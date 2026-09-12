% Permute CN->VN message -------------------------------------------------------
function C2V_p = permute_C2V(h, C2V)
global GF_MUL;
N_GF = 6;           % Number of GF(q) bits
Q_GF = 2^N_GF; 
    C2V_p = zeros(1, Q_GF, 'single');
    for i = 1:Q_GF
        C2V_p(i) = C2V(GF_MUL(h+1, i)+1);
    end
end