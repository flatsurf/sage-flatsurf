# Launch or Reinstall SageMath
$ErrorActionPreference = "Stop"

$appDataDir = [System.Environment]::GetFolderPath("LocalApplicationData")
$localDir = Join-Path -Path $appDataDir -ChildPath "sage-flatsurf-VERSION"
$wsldlExe = Join-Path -Path $localDir -ChildPath "sage-flatsurf-VERSION.exe"

Push-Location $appDataDir
$winHome = $(wsl wslpath -a ($HOME -replace '\\', '/')).Trim()
Pop-Location

function Usage {
  Write-Output "Usage $($MyInvocation.MyCommand.Path) --jupyterlab | --repl | --shell | --reinstall | --uninstall | --install"
}

function EnsureInstalled {
  $consoleHandle = (Get-Process -Id $PID).MainWindowHandle

  if ($consoleHandle -ne 0) {
    if (Test-Path $wsldlExe) {
      Write-Output "Not installing, wsldl executable already exists"
      return
    }

    Install
  } else {
    # Restart script in a visible terminal
    Start-Process powershell.exe -ArgumentList @("-File", "`"$PSCommandPath`"", $commandArgs[0]) -Wait
    exit
  }
}

function Install {
  New-Item -ItemType Directory -Path $localDir -Force > $null
  foreach ($fileName in @("sage-flatsurf-VERSION.exe", "preset.json")) {
    $sourceDir = $PSScriptRoot
    $sourceFile = Join-Path -Path $sourceDir -ChildPath $fileName
    $destFile = Join-Path -Path $localDir -ChildPath $fileName
    Copy-Item -Path $sourceFile -Destination $destFile -Force
  }

  Write-Output "Installing Linux base system"

  Push-Location $localDir
  # Ignore the request to press enter to continue.
  echo "" | & "$wsldlExe"

  Write-Output "We already pressed Enter for you."

  Write-Output "Setting up browser support"

  & "$wsldlExe" "run" "apt" "update"
  & "$wsldlExe" "run" "apt" "install" "-y" "wslu" "xdg-utils"

  Write-Output "Setting hostname"

  & "$wsldlExe" "run" "hostnamectl" "set-hostname" "sage-flatsurf-VERSION"

  Write-Output "Setting up wsl user"

  & "$wsldlExe" "run" "adduser" "wsl" "--disabled-password" "--gecos" "''"
  & "$wsldlExe" "config" "--default-user" "wsl"

  Write-Output "Preparing pixi installer"

  $tarball = Join-Path -Path $sourceDir -ChildPath "sage-flatsurf-VERSION.unix.tar.gz"
  Copy-Item -Path $tarball -Destination "\\wsl$\sage-flatsurf-VERSION\home\wsl\sage-flatsurf-VERSION.unix.tar.gz"

  & "$wsldlExe" "run" "sh" "-c" "cd ~ && tar zxf sage-flatsurf-VERSION.unix.tar.gz"
  Pop-Location
}

function JupyterLab {
  EnsureInstalled
  Start-Process -FilePath "$wsldlExe" -ArgumentList "run", "sh", "-c", "`"JUPYTERLAB_HOME=$winHome /home/wsl/sage-flatsurf-VERSION/jupyterlab`""
}

function REPL {
  EnsureInstalled
  Start-Process -FilePath "$wsldlExe" -ArgumentList "run", "/home/wsl/sage-flatsurf-VERSION/sage"
}

function Shell {
  EnsureInstalled
  Start-Process -FilePath "$wsldlExe" -ArgumentList "run", "/home/wsl/sage-flatsurf-VERSION/shell"
}

function Uninstall {
  try {
    & wsl --unregister sage-flatsurf-VERSION
  } catch {
    Write-Output "Ignoring failed unregister of VM, it is probably not registered yet."
  }
  
  Remove-Item -Path $localDir -Recurse -Force
}

function Reinstall {
  Uninstall
  Install
}

if ($args.Count -ne 1) {
  Usage
  exit 1
}

$commandArgs = $args

switch ($args[0]) {
  '--jupyterlab' {
    Write-Output "Launching JupyterLab..."
    JupyterLab
  }
  '--repl' {
    Write-Output "Launching SageMath REPL..."
    REPL
  }
  '--shell' {
    Write-Output "Launching a WSL shell..."
    Shell
  }
  '--install' {
    Write-Output "Installing virtual machine..."
    Install
  }
  '--reinstall' {
    Write-Output "Reinstalling virtual machine..."
    Reinstall
  }
  '--uninstall' {
    Write-Output "Removing virtual machine..."
    Uninstall
  }
  default {
    Usage
    exit 1
  }
}
