Add-Type -AssemblyName System.Windows.Forms

$result = [System.Windows.Forms.MessageBox]::Show("Your version of Windows Subsystem for Linux might not be up-to-date. The latest version is required. Do you want to update?", "Confirmation", [System.Windows.Forms.MessageBoxButtons]::YesNo)

if ($result -eq [System.Windows.Forms.DialogResult]::No) {
    exit
}

$kernelInstallerUrl = "https://wslstorestorage.blob.core.windows.net/wslblob/wsl_update_x64.msi"
$installerPath = "$env:TEMP\wsl_update_x64.msi"
Invoke-WebRequest -Uri $kernelInstallerUrl -OutFile $installerPath
Start-Process msiexec.exe -ArgumentList "/i", "`"$installerPath`"", "/quiet", "/norestart" -Wait
Remove-Item $installerPath -Force
Write-Host "WSL kernel has been updated successfully."
