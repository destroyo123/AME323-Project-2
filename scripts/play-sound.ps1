param(
    [string]$SoundFile
)

Add-Type -AssemblyName presentationCore
$player = New-Object System.Windows.Media.MediaPlayer
$player.Open([Uri]::new((Resolve-Path $SoundFile)))
$player.Play()

Start-Sleep -Seconds 3