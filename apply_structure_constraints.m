function L_constraints = apply_structure_constraints(L, N, m_sub, n_sub, ContStruc)
    % Based on Page 12: Imposing structure to L identical to Kx [cite: 259]
    L_constraints = [];
    row_ptr = 0;
    for i = 1:N
        col_ptr = 0;
        for j = 1:N
            % If the communication link (j -> i) does not exist in ContStruc
            if ContStruc(i,j) == 0
                % Constrain the corresponding block in L to be zero [cite: 242, 259]
                L_constraints = [L_constraints, ...
                    L(row_ptr+1:row_ptr+m_sub(i), col_ptr+1:col_ptr+n_sub(j)) == 0];
            end
            col_ptr = col_ptr + n_sub(j);
        end
        row_ptr = row_ptr + m_sub(i);
    end
end