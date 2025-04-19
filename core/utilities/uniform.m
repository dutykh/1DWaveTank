function [xc, dx, x_edge] = uniform(domain, N)
    
    %UNIFORM  Generate one‑dimensional uniform mesh.
    %   [xc,dx,x_edge] = UNIFORM(domain, N) returns the cell centres XC (column
    %   vector, length N), uniform mesh‑spacing DX and the vector of cell edges
    %   X_EDGE (length N+1) spanning [domain.xmin, domain.xmax].
    
        validateattributes(domain.xmin, {"numeric"}, {"scalar", "real"});
        validateattributes(domain.xmax, {"numeric"}, {"scalar", "real", ">", domain.xmin});
        validateattributes(N,            {"numeric"}, {"scalar", "integer", ">", 2});
    
        x_edge = linspace(domain.xmin, domain.xmax, N+1).';
        dx     = x_edge(2) - x_edge(1);
        xc     = 0.5*(x_edge(1:end-1) + x_edge(2:end));

end