param (
    [switch]$SkipTests = $false
)

$ErrorActionPreference = "Stop"

$repoRoot = $PSScriptRoot
if ([string]::IsNullOrEmpty($repoRoot)) {
    $repoRoot = Get-Location
}

Write-Host "========================================================" -ForegroundColor Cyan
Write-Host "Starting End-to-End Analysis Pipeline" -ForegroundColor Cyan
Write-Host "========================================================" -ForegroundColor Cyan

# 1. Run MATLAB data processing pipeline
Write-Host "`n[1/2] Running MATLAB workflows (execute_all_workflows.m)..." -ForegroundColor Yellow

$matlabCmd = "try, run('$repoRoot\execute_all_workflows.m'); catch ME, disp(ME.message); exit(1); end; exit(0);"
if ($SkipTests) {
    Write-Host "Skipping MATLAB tests as requested..." -ForegroundColor DarkGray
    # Temporarily modify config to skip tests if requested via parameter
    $configPath = Join-Path $repoRoot "config.json"
    if (Test-Path $configPath) {
        $config = Get-Content $configPath | ConvertFrom-Json
        $originalSkip = $null
        if ($config.PSObject.Properties.Match('skip_tests').Count -gt 0) {
            $originalSkip = $config.skip_tests
            $config.skip_tests = $true
        } else {
            $config | Add-Member -MemberType NoteProperty -Name "skip_tests" -Value $true
        }
        $config | ConvertTo-Json -Depth 10 | Set-Content $configPath
    }
}

try {
    # Check if MATLAB is in PATH
    $matlabPath = Get-Command "matlab" -ErrorAction SilentlyContinue
    if ($null -eq $matlabPath) {
        Write-Host "WARNING: 'matlab' is not in your system PATH." -ForegroundColor Red
        Write-Host "Please run 'execute_all_workflows.m' manually in MATLAB first," -ForegroundColor Red
        Write-Host "then run this script again." -ForegroundColor Red
        # Pause briefly then proceed to try the Python pipeline anyway on the latest existing folder
        Start-Sleep -Seconds 2
    } else {
        # Run MATLAB in batch mode without desktop
        Start-Process -FilePath "matlab" -ArgumentList "-batch `"$matlabCmd`"" -Wait -NoNewWindow
        if ($LASTEXITCODE -ne 0) {
            throw "MATLAB pipeline failed with exit code $LASTEXITCODE"
        }
        Write-Host "MATLAB pipeline completed successfully." -ForegroundColor Green
    }
} catch {
    Write-Host "Error launching MATLAB: $_" -ForegroundColor Red
} finally {
    if ($SkipTests -and (Test-Path $configPath)) {
        # Restore test setting
        $config = Get-Content $configPath | ConvertFrom-Json
        if ($null -ne $originalSkip) {
            $config.skip_tests = $originalSkip
        } else {
            $config.PSObject.Properties.Remove('skip_tests')
        }
        $config | ConvertTo-Json -Depth 10 | Set-Content $configPath
    }
}

# 2. Find the most recently created saved_files directory
$latestSavedDir = Get-ChildItem -Path $repoRoot -Filter "saved_files_*" -Directory | Sort-Object CreationTime -Descending | Select-Object -First 1

if ($null -eq $latestSavedDir) {
    throw "Could not find any 'saved_files_*' output directory."
}
$savedPath = $latestSavedDir.FullName
Write-Host "Detected latest output directory: $($latestSavedDir.Name)" -ForegroundColor DarkGray


# 3. Run Python analysis pipeline
Write-Host "`n[2/2] Running Python analysis and HTML report generation..." -ForegroundColor Yellow

$analysisDir = Join-Path $repoRoot "analysis"
Set-Location $analysisDir

try {
    # 3a. Generate Python graphs and statistics
    Write-Host "Generating graph plots..." -ForegroundColor DarkGray
    python run_analysis.py --folder $savedPath
    if ($LASTEXITCODE -ne 0) { throw "Graph generation failed." }

    # 3b. Parse MAT metrics files
    Write-Host "Parsing MAT metrics..." -ForegroundColor DarkGray
    python parse_mat_metrics.py $savedPath
    if ($LASTEXITCODE -ne 0) { throw "MAT metrics parsing failed." }

    # 3c. Generate final HTML report
    Write-Host "Building HTML report..." -ForegroundColor DarkGray
    python generate_report.py $savedPath
    if ($LASTEXITCODE -ne 0) { throw "HTML report generation failed." }
} catch {
    Set-Location $repoRoot
    throw $_
}

Set-Location $repoRoot


Write-Host "`n========================================================" -ForegroundColor Cyan
Write-Host "Pipeline Complete!" -ForegroundColor Green
Write-Host "Review your report at: $($latestSavedDir.FullName)\analysis_report.html" -ForegroundColor Cyan
Write-Host "========================================================" -ForegroundColor Cyan
