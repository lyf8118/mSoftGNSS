function V2C_p = permute_V2C(h, V2C)
global GF_MUL;
N_GF = 6;           % Number of GF(q) bits
Q_GF = 2^N_GF; 
    V2C_p = zeros(1, Q_GF, 'single');
    for i = 1:Q_GF
         % 将 MATLAB 索引映射到 GF 元素（0-based）
        V2C_p(GF_MUL(h+1, i)+1) = V2C(i);
    end
end