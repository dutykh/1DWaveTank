function F = PhysFlux(v, g)

    %PHYSFLUX   Computes the physical flux for the shallow water equations.
    %   F = PhysFlux(v, g) returns the physical flux for state v and gravity g.
    
    h = v(:,1);
    hu = v(:,2);
    u = hu ./ (h + eps);
    F = [hu, hu .* u + 0.5 * g * h.^2];

end