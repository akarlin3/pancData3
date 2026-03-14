function write_sentinel_file(output_folder, prefix, message, dwi_type_name)
% WRITE_SENTINEL_FILE — Write a sentinel file confirming step completion.
%
%   Creates a text file at
%     output_folder/<prefix>_<dwi_type_name>.txt
%   containing the given message.  Used by analysis scripts and test
%   harnesses to verify that a pipeline step ran successfully.
%
%   Parameters
%   ----------
%   output_folder : char
%       Directory in which to write the sentinel file.
%   prefix : char
%       File name prefix (e.g. 'metrics_longitudinal_results').
%   message : char
%       Message to write into the file (e.g. 'Longitudinal metrics generated successfully.').
%   dwi_type_name : char
%       DWI type label for file naming (e.g. 'Standard').
%
%   See also: run_dwi_pipeline, execute_pipeline_step

    sentinel_file = fullfile(output_folder, sprintf('%s_%s.txt', prefix, dwi_type_name));
    fid = fopen(sentinel_file, 'w');
    if fid < 0
        warning('write_sentinel_file:fileWriteFailed', 'Cannot write %s', sentinel_file);
    else
        fprintf(fid, '%s\n', message);
        fclose(fid);
    end
    fprintf('      💾 Saved %s to %s\n', prefix, sentinel_file);
end
