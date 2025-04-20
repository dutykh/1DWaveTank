function status=odetpbar(t,y,flag)

    persistent tf tstart;
    
    if isempty(flag)
        % Integration steps
        ts=mean(t);
        progress=100*ts/tf;
        core.utils.textprogressbar(progress); % Use namespaced call
        status=0;
    else
        switch flag
            case 'init'     % Initializing progress bar
                tstart=tic;
                tf=max(t);
                core.utils.textprogressbar('ODE integration: '); % Use namespaced call
            case 'done'     % Finishing status function
                tf=[];
                core.utils.textprogressbar(''); % Use namespaced call
                fprintf('\n   Integration time: %.3f s\n', toc(tstart)); % Added newline and formatting
                tstart=[];
            otherwise
                error('odetpbar:UnknownError',...
                    'Unknown error has occured');
        end
        status = 0; % Ensure status is always returned
    end

end