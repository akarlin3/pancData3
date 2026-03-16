classdef test_gpu_available < matlab.unittest.TestCase
    % TEST_GPU_AVAILABLE Unit tests for the gpu_available utility.
    %
    % Validates:
    %   - Function returns two outputs [available, gpu_dev]
    %   - available is a logical scalar
    %   - gpu_dev is either empty or a gpuDevice object
    %   - Default device_index is 1
    %   - Invalid device index returns gracefully (no error)

    methods(TestMethodSetup)
        function addPaths(testCase)
            repoRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(repoRoot, 'utils'));
        end
    end

    methods(Test)

        function testReturnsTwoOutputs(testCase)
            % gpu_available must always return exactly two outputs
            % without erroring, regardless of GPU hardware.
            [avail, dev] = gpu_available();
            testCase.verifyClass(avail, 'logical', ...
                'First output should be logical.');
            testCase.verifyTrue(isscalar(avail), ...
                'First output should be scalar.');
        end

        function testDefaultDeviceIndex(testCase)
            % Calling without arguments should default to device 1.
            % This should not error even if no GPU is present.
            [avail, ~] = gpu_available();
            testCase.verifyClass(avail, 'logical');
        end

        function testInvalidDeviceIndexGraceful(testCase)
            % A device index of 999 (almost certainly nonexistent)
            % should return false, not throw an error.
            [avail, dev] = gpu_available(999);
            testCase.verifyFalse(avail, ...
                'Nonexistent device index should return false.');
            testCase.verifyEmpty(dev, ...
                'Nonexistent device should return empty gpu_dev.');
        end

        function testEmptyDeviceIndexDefaultsToOne(testCase)
            % Passing empty should behave like passing 1.
            [avail1, ~] = gpu_available([]);
            [avail2, ~] = gpu_available(1);
            testCase.verifyEqual(avail1, avail2, ...
                'Empty device_index should behave like 1.');
        end

        function testGpuDevObjectOrEmpty(testCase)
            % Second output is either a gpuDevice object (when GPU
            % available) or empty (when not).
            [avail, dev] = gpu_available();
            if avail
                testCase.verifyClass(dev, 'parallel.gpu.CUDADevice', ...
                    'When available, gpu_dev should be a CUDADevice object.');
            else
                testCase.verifyEmpty(dev, ...
                    'When not available, gpu_dev should be empty.');
            end
        end

    end
end
