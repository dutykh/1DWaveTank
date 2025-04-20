% +core/+utils/get_bc_style.m
function style = get_bc_style(handle)

    handle_name = func2str(handle);
    if contains(handle_name, 'wall')
        style.color = 'k'; style.style = '-'; % Solid black for wall
    elseif contains(handle_name, 'generating')
        style.color = 'b'; style.style = '--'; % Dashed blue for generating
    elseif contains(handle_name, 'open')
        style.color = [0.5 0.5 0.5]; style.style = ':'; % Dotted gray for open
    else
        style.color = 'm'; style.style = '-.'; % Magenta dash-dot for unknown/other
    end

end