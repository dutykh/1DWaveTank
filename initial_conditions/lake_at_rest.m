function U0 = lake_at_rest(x, cfg)
    
    %LAKE_AT_REST  Zero‑surface‑elevation, zero‑velocity initial state.
    %   For NSWE the conserved variables are [h; hu].  Here h = H0 and hu = 0.
    
        h  = cfg.param.H0 * ones(size(x));
        hu = zeros(size(x));
        U0 = [h.'; hu.'];            % row‑major to comply with solver conventions

end