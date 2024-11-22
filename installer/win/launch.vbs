Set objArgs = WScript.Arguments
Set objFSO = CreateObject("Scripting.FileSystemObject")
scriptDir = objFSO.GetParentFolderName(WScript.ScriptFullName)
psScriptPath = scriptDir & "\launch.ps1"

Dim psArgs, windowStyle, windowPolicy, waitPolicy
psArgs = ""
windowStyle = "Normal"
windowPolicy = 1
waitPolicy = True

For i = 0 To objArgs.Count - 1
    If LCase(objArgs(i)) = "--quiet" Then
        windowStyle = "Hidden"
        windowPolicy = 0
    Else
        psArgs = psArgs & " """ & objArgs(i) & """"
    End If
Next

Set objShell = CreateObject("WScript.Shell")
objShell.Run "powershell.exe -NoProfile -ExecutionPolicy Bypass -WindowStyle " & windowStyle & " -File """ & psScriptPath & """" & " --ensure-install", windowPolicy, waitPolicy
Set objShell = Nothing

Set objShell = CreateObject("WScript.Shell")
objShell.Run "powershell.exe -NoProfile -ExecutionPolicy Bypass -WindowStyle " & windowStyle & " -File """ & psScriptPath & """" & psArgs, windowPolicy, waitPolicy
Set objShell = Nothing

Set objFSO = Nothing
Set objShell = Nothing
