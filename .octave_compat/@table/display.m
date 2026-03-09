% DISPLAY  Console display method for the Octave table shim.
%
%   Replaces MATLAB's table display formatting. MATLAB's table prints a
%   nicely formatted columnar layout with headers; this shim simply dumps
%   the underlying struct fields since Octave lacks the rich table renderer.
function display(obj)
    % Print a label and then delegate to struct display for simplicity.
    disp('table object:');
    disp(struct(obj));
end
