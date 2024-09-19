[Setup]
AppName=flatsurf
AppVersion=0.5.2
WizardStyle=modern
DefaultDirName={autopf}\sage-flatsurf-0.5.2
DefaultGroupName=sage-flatsurf-0.5.2
OutputBaseFilename=sage-flatsurf-0.5.2
SetupIconFile=sage-flatsurf-0.5.2.ico
UninstallDisplayIcon=sage-flatsurf-0.5.2.ico
; TODO
;UninstallDisplayIcon={app}\MyProg.exe

; TODO autostart after install, see https://stackoverflow.com/q/3767227/812379

; TODO: Have a look at all the [Setup] section directives in inno setup documentation.

[Files]
Source: "sage-flatsurf-0.5.2.exe"; DestDir: "{app}";
Source: "sage-flatsurf-0.5.2.ico"; DestDir: "{app}";
Source: "sage-flatsurf-0.5.2.tar.gz"; DestDir: "{app}";
Source: "preset.json"; DestDir: "{app}";
Source: "launch.ps1"; DestDir: "{app}";
Source: "launch.vbs"; DestDir: "{app}";
Source: "UpdateWSLKernel.ps1"; DestDir: "{app}";

[Icons]
Name: "{commondesktop}\SageMath (sage-flatsurf 0.5.2)"; Filename: "wscript.exe"; Parameters: """{app}\launch.vbs"" --repl --quiet"; WorkingDir: "{app}"; IconFilename: "{app}\sage-flatsurf-0.5.2.ico"; Tasks: desktopicon
Name: "{commondesktop}\JupyterLab (sage-flatsurf 0.5.2)"; Filename: "wscript.exe"; Parameters: """{app}\launch.vbs"" --jupyterlab --quiet"; WorkingDir: "{app}"; IconFilename: "{app}\sage-flatsurf-0.5.2.ico"; Tasks: desktopicon
Name: "{commondesktop}\Shell (sage-flatsurf 0.5.2)"; Filename: "wscript.exe"; Parameters: """{app}\launch.vbs"" --shell --quiet"; WorkingDir: "{app}"; IconFilename: "{app}\sage-flatsurf-0.5.2.ico"; Tasks: desktopicon
Name: "{commonprograms}\sage-flatsurf 0.5.2\SageMath"; Filename: "wscript.exe"; Parameters: """{app}\launch.vbs"" --repl --quiet"; WorkingDir: "{app}"; IconFilename: "{app}\sage-flatsurf-0.5.2.ico"
Name: "{commonprograms}\sage-flatsurf 0.5.2\JupyterLab"; Filename: "wscript.exe"; Parameters: """{app}\launch.vbs"" --jupyterlab --quiet"; WorkingDir: "{app}"; IconFilename: "{app}\sage-flatsurf-0.5.2.ico"
Name: "{commonprograms}\sage-flatsurf 0.5.2\Shell"; Filename: "wscript.exe"; Parameters: """{app}\launch.vbs"" --shell --quiet"; WorkingDir: "{app}"; IconFilename: "{app}\sage-flatsurf-0.5.2.ico"
Name: "{commonprograms}\sage-flatsurf 0.5.2\Uninstall Virtual Machine"; Filename: "wscript.exe"; Parameters: """{app}\launch.vbs"" --uninstall"; WorkingDir: "{app}"; IconFilename: "{app}\sage-flatsurf-0.5.2.ico"
Name: "{commonprograms}\sage-flatsurf 0.5.2\Reinstall Virtual Machine"; Filename: "wscript.exe"; Parameters: """{app}\launch.vbs"" --reinstall"; WorkingDir: "{app}"; IconFilename: "{app}\sage-flatsurf-0.5.2.ico"


[Tasks]
Name: "desktopicon"; Description: "Create &desktop shortcuts"; GroupDescription: "Additional icons:"; Flags: unchecked

[Run]
Filename: "powershell.exe"; Description: "Update Linux Kernel (often required)"; Parameters: "-ExecutionPolicy Bypass -NoProfile -File ""{app}\UpdateWSLKernel.ps1"""; Flags: postinstall

[Code]
const
  RunOnceName = 'sage-flatsurf-0.5.2 restart';

  QuitMessageError = 'Error. Cannot continue.';

var
  Restarted: Boolean;

// Complain if the user didn't restart as requested, otherwise continue.
function InitializeSetup(): Boolean;
begin
  Restarted := ExpandConstant('{param:restart|0}') = '1';

  if not Restarted then begin
    Result := not RegValueExists(HKA, 'Software\Microsoft\Windows\CurrentVersion\RunOnce', RunOnceName);
    if not Result then
      MsgBox('Please restart your machine to finish the installation.', mbError, mb_Ok);
  end else
    Result := True;
end;

// Read the string from redirect and trim trailing whitespace
function ReadStdout(const redirect: String): String;
var
  stdout_ansi : AnsiString;
begin
  LoadStringFromFile(redirect, stdout_ansi);
  DeleteFile(redirect);
  Result := String(stdout_ansi);
  if (Length(Result) >= 2) and (Result[Length(Result) - 1] = #13) and (Result[Length(Result)] = #10) then
    Delete(Result, Length(Result) - 1, 2);
end;

// Return whether WSL2 is installed, return "Enabled", "Disabled", or an error message.
function DetectWSL: String;
var
  ResultCode: Integer;
  redirect : String;
  stdout : String;
begin
  redirect := ExpandConstant('{tmp}\~sage-flatsurf.stdout');

  if not Exec('powershell', '-NoProfile -ExecutionPolicy Bypass -Command "(Get-WindowsOptionalFeature -Online -FeatureName Microsoft-Windows-Subsystem-Linux).State | out-file -encoding ASCII \"' + redirect + '\""', '', SW_HIDE, ewWaitUntilTerminated, ResultCode) then begin
    Result := 'Could not determine whether Windows Subsystem for Linux is already installed. The installation will now abort.';
    Exit;
  end;
  stdout := ReadStdout(redirect);
  Log('WSL State = ' + stdout);

  // With French culture set, this also said Enabled and Disabled, so apparently these are not localized.
  if stdout = 'Disabled' then begin
    Result := 'Disabled';
    Exit;
  end;
  if stdout <> 'Enabled' then begin
    Result := 'Could not determine whether Windows Subsystem for Linux is already installed: ' + stdout + '. The installation will now abort.';
    Exit;
  end;

  if not Exec('powershell', '-NoProfile -ExecutionPolicy Bypass -Command "(Get-WindowsOptionalFeature -Online -FeatureName VirtualMachinePlatform).State | out-file -encoding ASCII \"' + redirect + '\""', '', SW_HIDE, ewWaitUntilTerminated, ResultCode) then begin
    Result := 'Could not determine whether the Virtual Machine Platform is already installed. The installation will now abort.';
    Exit;
  end;
  stdout := ReadStdout(redirect);
  Log('WSL2 State = ' + stdout);

  // With French culture set, this also said Enabled and Disabled, so apparently these are not localized.
  if stdout = 'Disabled' then begin
    Result := 'Disabled';
    Exit;
  end;
  if stdout <> 'Enabled' then begin
    Result := 'Could not determine whether the Virtual Machine Platform is already installed: ' + stdout + '. The installation will now abort.';
    Exit;
  end;

  Result := 'Enabled';
end;

// Return whether WSL2 could be installed, return "Enabled" or an error message.
function EnableWSL: String;
var
  ResultCode: Integer;
begin
  Log('Enabling WSL')
  if not Exec('powershell', '-NoProfile -ExecutionPolicy Bypass -Command "Enable-WindowsOptionalFeature -Online -FeatureName Microsoft-Windows-Subsystem-Linux -NoRestart"', '', SW_HIDE, ewWaitUntilTerminated, ResultCode) then begin
    Result := 'Could not enable Windows Subsystem for Linux. Cannot proceed installation without WSL.';
    Exit;
  end;

  if not Exec('powershell', '-NoProfile -ExecutionPolicy Bypass -Command "Enable-WindowsOptionalFeature -Online -FeatureName VirtualMachinePlatform -NoRestart"', '', SW_HIDE, ewWaitUntilTerminated, ResultCode) then begin
    Result := 'Could not enable the Virtual Machine Platform. Cannot proceed installation without WSL2.';
    Exit;
  end;

  Result := 'Enabled';
end;

function Quote(const S: String): String;
begin
  Result := '"' + S + '"';
end;

function AddParam(const S, P, V: String): String;
begin
  if V <> '""' then
    Result := S + ' /' + P + '=' + V;
end;

function AddSimpleParam(const S, P: String): String;
begin
 Result := S + ' /' + P;
end;

procedure CreateRunOnceEntry;
var
  RunOnceData: String;
begin
  RunOnceData := Quote(ExpandConstant('{srcexe}')) + ' /restart=1';
  RunOnceData := AddParam(RunOnceData, 'LANG', ExpandConstant('{language}'));
  RunOnceData := AddParam(RunOnceData, 'DIR', Quote(WizardDirValue));
  RunOnceData := AddParam(RunOnceData, 'GROUP', Quote(WizardGroupValue));
  if WizardNoIcons then
    RunOnceData := AddSimpleParam(RunOnceData, 'NOICONS');
  RunOnceData := AddParam(RunOnceData, 'TYPE', Quote(WizardSetupType(False)));
  RunOnceData := AddParam(RunOnceData, 'COMPONENTS', Quote(WizardSelectedComponents(False)));
  RunOnceData := AddParam(RunOnceData, 'TASKS', Quote(WizardSelectedTasks(False)));

  RegWriteStringValue(HKA, 'Software\Microsoft\Windows\CurrentVersion\RunOnce', RunOnceName, RunOnceData);
end;

function PrepareToInstall(var NeedsRestart: Boolean): String;
begin
  Result := DetectWSL;
  Log('DetectWSL = ' + Result);
  if Result = 'Disabled' then begin
    // To get into this code path again run either of:
    // Disable-WindowsOptionalFeature -Online -FeatureName Microsoft-Windows-Subsystem-Linux -NoRestart
    // Disable-WindowsOptionalFeature -Online -FeatureName VirtualMachinePlatform -NoRestart  
    if MsgBox('Windows Subsystem for Linux is required but not configured (correctly) on your system. Do you want me to reconfigure it?', mbConfirmation, MB_YESNO) = IDNO then begin
      Result := 'Cannot finish installation without Windows Subsystem for Linux 2';
      Exit;
    end;
    EnableWSL;
    CreateRunOnceEntry;
    NeedsRestart := True;
    Result := 'A reboot is required to finish the installation.';
  end;

  if Result = 'Enabled' then begin
    Result := '';
  end;
end;

function ShouldSkipPage(PageID: Integer): Boolean;
begin
  Result := Restarted;
end;

procedure CurUninstallStepChanged(CurUninstallStep: TUninstallStep);
begin
  if CurUninstallStep = usUninstall then
  begin
    MsgBox('Make sure to run the "Uninstall Virtual Machine" entry from the Start menu before uninstalling sage-flatsurf. Otherwise, there might be a leftover virtual machine in your %LOCALAPPDATA% directory (which you can delete safely.)', mbConfirmation, MB_OK);
  end;
end;