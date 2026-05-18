Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

Set-Location -LiteralPath $PSScriptRoot

if (-not (Get-Command Rscript -ErrorAction SilentlyContinue)) {
    throw "Rscript was not found on PATH. Install R or add Rscript.exe to PATH, then rerun this script."
}

Rscript run_all.R

