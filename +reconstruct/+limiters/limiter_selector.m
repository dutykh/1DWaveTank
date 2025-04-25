function limiter_handle = limiter_selector(limiter_name)
% limiter_selector  Return a function handle to the requested MUSCL limiter
%
% Usage:
%   handle = reconstruct.limiters.limiter_selector('vanleer');
%   handle = reconstruct.limiters.limiter_selector('minmod');
%   ...
%
% Supported limiters: 'minmod', 'vanleer', 'superbee', 'mc', 'koren', 'umist', 'ospre', 'sweby', 'vanalbada'
% (case-insensitive)
%
% Returns @reconstruct.limiters.minmod, etc., or throws an error if unknown.

    if nargin < 1 || isempty(limiter_name)
        limiter_name = 'minmod';
    end
    name = lower(limiter_name);
    switch name
        case 'minmod'
            limiter_handle = @reconstruct.limiters.minmod;
        case 'vanleer'
            limiter_handle = @reconstruct.limiters.vanleer;
        case 'superbee'
            limiter_handle = @reconstruct.limiters.superbee;
        case 'mc'
            limiter_handle = @reconstruct.limiters.mc;
        case 'koren'
            limiter_handle = @reconstruct.limiters.koren;
        case 'umist'
            limiter_handle = @reconstruct.limiters.umist;
        case 'ospre'
            limiter_handle = @reconstruct.limiters.ospre;
        case 'sweby'
            limiter_handle = @reconstruct.limiters.sweby;
        case 'vanalbada'
            limiter_handle = @reconstruct.limiters.vanalbada;
        otherwise
            error('Unknown limiter "%s". Supported: minmod, vanleer, superbee, mc, koren, umist, ospre, sweby, vanalbada', limiter_name);
    end
end
