param(
    [string]$SoundFile
)

$player = New-Object System.Media.SoundPlayer $SoundFile
$player.PlaySync()