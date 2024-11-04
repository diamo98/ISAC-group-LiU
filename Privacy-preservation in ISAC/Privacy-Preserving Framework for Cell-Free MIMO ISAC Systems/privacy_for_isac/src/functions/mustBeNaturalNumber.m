function mustBeNaturalNumber(value)
    % check if finite natural number (including zero)
    if value < 0 || mod(value,1) ~= 0 || abs(value) == inf
        error("Value must be natural number")
    end
end

