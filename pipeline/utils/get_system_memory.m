function [total_gb, avail_gb] = get_system_memory()
% GET_SYSTEM_MEMORY  Return total and available system memory in GB.
%
%   [total_gb, avail_gb] = get_system_memory()
%
%   Uses platform-specific commands to query physical memory.
%   Returns [NaN, NaN] if memory cannot be determined (e.g., sandboxed
%   environments, unsupported platforms, or command failures).
%
%   See also: load_dwi_data, memory

    total_gb = NaN;
    avail_gb = NaN;

    try
        if ispc
            % Windows: use built-in memory() function
            [usr, sys] = memory;
            total_gb = sys.PhysicalMemory.Total / 1e9;
            avail_gb = sys.PhysicalMemory.Available / 1e9;
        elseif ismac
            % macOS: use sysctl
            [~, out] = system('sysctl -n hw.memsize');
            total_bytes = str2double(strtrim(out));
            if ~isnan(total_bytes)
                total_gb = total_bytes / 1e9;
            end
            % Available memory: use vm_stat
            [~, vm_out] = system('vm_stat');
            page_size = 4096;
            free_match = regexp(vm_out, 'Pages free:\s+(\d+)', 'tokens');
            inactive_match = regexp(vm_out, 'Pages inactive:\s+(\d+)', 'tokens');
            free_pages = 0;
            if ~isempty(free_match), free_pages = free_pages + str2double(free_match{1}{1}); end
            if ~isempty(inactive_match), free_pages = free_pages + str2double(inactive_match{1}{1}); end
            avail_gb = free_pages * page_size / 1e9;
        else
            % Linux: use /proc/meminfo
            [~, out] = system('cat /proc/meminfo');
            total_match = regexp(out, 'MemTotal:\s+(\d+)', 'tokens');
            avail_match = regexp(out, 'MemAvailable:\s+(\d+)', 'tokens');
            if ~isempty(total_match)
                total_gb = str2double(total_match{1}{1}) / 1e6; % KB to GB
            end
            if ~isempty(avail_match)
                avail_gb = str2double(avail_match{1}{1}) / 1e6;
            end
        end
    catch
        % Silently return NaN if memory query fails
    end
end
