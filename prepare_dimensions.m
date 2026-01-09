function [Btot, ntot, mtot, m_sub, n_sub] = prepare_dimensions(Bdec, Cdec, N)
    Btot = [];
    m_sub = zeros(1,N); n_sub = zeros(1,N);
    for i=1:N
        m_sub(i) = size(Bdec{i}, 2);
        n_sub(i) = size(Cdec{i}, 1);
        Btot = [Btot, Bdec{i}];
    end
    ntot = sum(n_sub);
    mtot = sum(m_sub);
end