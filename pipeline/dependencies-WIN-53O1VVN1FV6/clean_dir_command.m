% perform the "dir" command, but remove extraneous entries (ignore ., ..,
% and non-folder)
% EA 04.03.2025

function dir_cleaned = clean_dir_command(folder)
    dir_cleaned = dir(folder);
    dir_cleaned = dir_cleaned([dir_cleaned.isdir]);
    dir_cleaned = dir_cleaned(~strcmpi({dir_cleaned.name},'.'));
    dir_cleaned = dir_cleaned(~strcmpi({dir_cleaned.name},'..'));
end