classdef test_parse_clinical_value < matlab.unittest.TestCase
    % TEST_PARSE_CLINICAL_VALUE  Coercion of clinical-spreadsheet cells.
    %
    % Locks in behaviour for the type variants readtable can return for a
    % binary outcome column. Regression: the live MSK cohort had ~7 LF
    % events in the spreadsheet, but the pipeline read every patient as
    % lf=0 because readtable inferred a non-numeric column type. The
    % helper coerces all variants to a numeric scalar so the downstream
    % assignment to a numeric array works regardless.

    methods(TestMethodSetup)
        function addUtilsToPath(testCase) %#ok<INUSD>
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods(Test)
        function testNumericZero(testCase)
            testCase.verifyEqual(parse_clinical_value(0), 0);
        end

        function testNumericOne(testCase)
            testCase.verifyEqual(parse_clinical_value(1), 1);
        end

        function testLogicalTrueAndFalse(testCase)
            testCase.verifyEqual(parse_clinical_value(true), 1);
            testCase.verifyEqual(parse_clinical_value(false), 0);
        end

        function testInt32(testCase)
            testCase.verifyEqual(parse_clinical_value(int32(1)), 1);
        end

        function testCharOne(testCase)
            % readtable returns char arrays when a column has mixed text
            testCase.verifyEqual(parse_clinical_value('1'), 1);
        end

        function testCharZeroWithWhitespace(testCase)
            testCase.verifyEqual(parse_clinical_value(' 0 '), 0);
        end

        function testCharEmptyIsNaN(testCase)
            testCase.verifyTrue(isnan(parse_clinical_value('')));
        end

        function testStringScalarOne(testCase)
            testCase.verifyEqual(parse_clinical_value("1"), 1);
        end

        function testStringMissingIsNaN(testCase)
            testCase.verifyTrue(isnan(parse_clinical_value(string(missing))));
        end

        function testCellOfCharOne(testCase)
            % Most likely readtable variant: cell-of-char column produces
            % a 1x1 cell when sliced by integer index.
            testCase.verifyEqual(parse_clinical_value({'1'}), 1);
        end

        function testCellOfCharZero(testCase)
            testCase.verifyEqual(parse_clinical_value({'0'}), 0);
        end

        function testCellOfNumericOne(testCase)
            testCase.verifyEqual(parse_clinical_value({1}), 1);
        end

        function testEmptyCellIsNaN(testCase)
            testCase.verifyTrue(isnan(parse_clinical_value({})));
        end

        function testEmptyArrayIsNaN(testCase)
            testCase.verifyTrue(isnan(parse_clinical_value([])));
        end

        function testCategoricalOne(testCase)
            c = categorical("1");
            testCase.verifyEqual(parse_clinical_value(c), 1);
        end

        function testCategoricalMissingIsNaN(testCase)
            c = categorical(missing);
            testCase.verifyTrue(isnan(parse_clinical_value(c)));
        end

        function testNoArgsIsNaN(testCase)
            testCase.verifyTrue(isnan(parse_clinical_value()));
        end
    end
end
