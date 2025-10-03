function [mat_out] = vectomat(vec_in,M)
    %mat_out = zeros(split,3)
    mat_out = [vec_in(1:M), vec_in(M+1:2*M) , vec_in(2*M+1:end) ];
end
