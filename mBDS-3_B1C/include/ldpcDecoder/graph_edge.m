function [ie, je, he, ne] = graph_edge(H_idx, H_ele)
    % Count the total number of nonzero elements.
    total_len = sum(cellfun(@length, H_idx));

    % Pre-allocate edge arrays.
    ie = zeros(1, total_len);
    je = zeros(1, total_len);
    he = zeros(1, total_len);

    % Fill edge lists.
    idx = 1;
    for i = 1:length(H_idx)
        len = length(H_idx{i});
        rng = idx:(idx+len-1);

        ie(rng) = i;
        je(rng) = H_idx{i};
        he(rng) = H_ele{i};

        idx = idx + len;
    end

    ne = total_len;
end
