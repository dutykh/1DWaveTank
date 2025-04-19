function h = flat(x, cfg)
    
    %FLAT  Constant-depth bathymetry.
    %   h = FLAT(x, cfg) returns a vector of still-water depths evaluated at the
    %   cell centres X.  The constant depth is taken from cfg.param.H0 (default
    %   0.5 m if unspecified).
    
        if ~isfield(cfg.param, "H0"),  cfg.param.H0 = 0.5;  end
        h = cfg.param.H0 * ones(size(x));

end