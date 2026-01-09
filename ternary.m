function res = ternary(condition, true_val, false_val)
    % Basit üçlü operatör (MATLAB'da yerleşik değildir)
    if condition
        res = true_val;
    else
        res = false_val;
    end
end