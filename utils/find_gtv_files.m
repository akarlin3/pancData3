function [gtvp_path, gtvn_path] = find_gtv_files(fxfolder, dwii, pat_name)
    % FIND_GTV_FILES Helper function to locate GTV mask files
    %   Handles multiple naming conventions and repeat indices

    gtvp_path = '';
    gtvn_path = '';

    if ~contains(pat_name, 'two')
        % find associated GTV (need to avoid using the date in the filenames to query)
        gtvp_path = discover_gtv_file(fxfolder, '*GTV*', dwii);
    else
        % --- Special logic for patients with both GTVp and GTVn ---
        % Search multiple naming conventions for the primary pancreatic GTV (GTVp)
        gtvp_patterns = {'*GTV_MR', '*GTVp', '*GTV_panc*'};
        gtvp_path = discover_gtv_file(fxfolder, gtvp_patterns, dwii);

        % Search for nodal GTV masks: GTV_LN, GTVn, GTV_node
        gtvn_patterns = {'*GTV*LN', '*GTVn', '*GTV_node*'};
        gtvn_path = discover_gtv_file(fxfolder, gtvn_patterns, dwii);
    end
end
