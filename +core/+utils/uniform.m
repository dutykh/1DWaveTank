function [xc, dx, x_edge] = uniform(domain, N)

    %UNIFORM Creates a uniform 1D mesh.
    %   [xc, dx, x_edge] = UNIFORM(domain, N) generates a uniform mesh
    %   for the spatial domain defined by the structure 'domain' (with fields
    %   'xmin' and 'xmax') using N control volumes (cells).
    %
    %   Outputs:
    %       xc      - Vector (N x 1) of cell centre coordinates.
    %       dx      - Scalar spatial step size.
    %       x_edge  - Vector ((N+1) x 1) of cell edge coordinates.
    
    if ~isfield(domain, 'xmin') || ~isfield(domain, 'xmax')
        error('Domain structure must contain xmin and xmax fields.');
    end
    
    x_edge = linspace(domain.xmin, domain.xmax, N+1)'; % Cell edges
    dx = (domain.xmax - domain.xmin) / N;             % Spatial step
    xc = x_edge(1:N) + 0.5 * dx;                      % Cell centres
    
end