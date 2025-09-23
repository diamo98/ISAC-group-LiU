function [s, bits] = randQPSK(L)
% RANDQPSK  Generate L random Gray-coded QPSK symbols (unit average power).
% Returns:
%   s    : Lx1 complex QPSK symbols in { (±1 ± j)/√2 }
%   bits : Lx2 bits [bI bQ] used for mapping

    bits = randi([0 1], L, 2);          % random bits [bI bQ]
    I = 1 - 2*bits(:,1);                % map 0->+1, 1->-1
    Q = 1 - 2*bits(:,2);
    s = (I + 1j*Q) / sqrt(2);           % Gray-coded QPSK constellation
end