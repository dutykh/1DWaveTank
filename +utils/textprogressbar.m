%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +core/+utils/textprogressbar.m
%
% Purpose:
%   Creates and updates a text-based progress bar in the MATLAB command window.
%   Allows displaying the percentage completion of a task (e.g., ODE integration).
%   Handles initialization, percentage updates, and termination/cleanup.
%
% Syntax:
%   textprogressbar(c)
%
% Inputs:
%   c   - Can be:
%         [char] Text string: Used to initialize the progress bar (e.g., 'Processing:')
%                or to terminate it (e.g., ' Done.'). An empty string clears the bar.
%         [numeric] Percentage value (0-100): Updates the bar to show the current progress.
%
% Outputs:
%   (none) - Prints directly to the command window.
%
% Dependencies:
%   None (self-contained utility).
%
% References:
%   - Inspired by: http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/
%   - Original Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
%   - Modified/Integrated by: Dr. Denys Dutykh (Khalifa University of Science and Technology, Abu Dhabi)
% Date:   21 April 2025 (Modification Date)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function textprogressbar(c)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialization & Persistent State                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Persistent variable to store the backspace characters needed to
    % overwrite the previous progress bar state in the command window.
    persistent strCR;

    % --- Visualization Parameters ---
    strPercentageLength = 10;   % [integer] Fixed width for the percentage string (e.g., ' 50%     ')
    strDotsMaximum      = 10;   % [integer] Total number of dots in the progress bar visualization

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main Logic: Handle Initialization, Update, Termination     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % --- Case 1: Initialization Error ---
    if isempty(strCR) && ~ischar(c)
        % If the progress bar hasn't been initialized (strCR is empty)
        % and the input is not a string, it's an error.
        error('textprogressbar:InitializationError', 'The text progress bar must be initialized with a string first.');

    % --- Case 2: Initialization ---
    elseif isempty(strCR) && ischar(c)
        % If strCR is empty and input is a string, this is the initialization call.
        fprintf('%s', c); % Print the initialization string (e.g., 'Processing: ')
        strCR = -1;       % Set strCR to -1 to indicate the first update needs no backspaces.

    % --- Case 3: Termination ---
    elseif ~isempty(strCR) && ischar(c)
        % If strCR is not empty (initialized) and input is a string, this is termination.
        strCR = []; % Clear the persistent variable to reset the state.
        fprintf([c '\n']); % Print the termination string (e.g., ' Done.') and a newline.

    % --- Case 4: Progress Update ---
    elseif isnumeric(c)
        % If input is numeric, update the progress bar.
        c = floor(c); % Ensure percentage is an integer.
        
        % Format the percentage string (fixed width)
        percentageOut = [num2str(c) '%%'];
        percentageOut = [percentageOut repmat(' ', 1, strPercentageLength - length(percentageOut) - 1)];
        
        % Format the dot bar string
        nDots = floor(c / 100 * strDotsMaximum); % Calculate number of dots to display
        dotOut = ['[' repmat('.', 1, nDots) repmat(' ', 1, strDotsMaximum - nDots) ']'];
        
        % Combine parts into the final string
        strOut = [percentageOut dotOut];

        % --- Print to Command Window ---
        % Use backspaces (stored in strCR) to overwrite the previous bar.
        if strCR == -1
            % First update: Just print, no backspaces needed.
            fprintf(strOut);
        else
            % Subsequent updates: Print backspaces then the new bar.
            fprintf([strCR strOut]);
        end

        % --- Update Backspace String ---
        % Store the correct number of backspace characters ('\b') needed
        % to erase the currently printed bar (strOut) on the next update.
        strCR = repmat('\b', 1, length(strOut) - 1);

    % --- Case 5: Error ---
    else
        % Handle any other unexpected input type.
        error('textprogressbar:InputError', 'Unsupported argument type');
    end

end
