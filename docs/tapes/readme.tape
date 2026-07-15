# VHS tape: end-to-end showcase for the README / docs landing page.
# One continuous session (no `clear`) so it reads like a real terminal.
Output docs/images/readme.gif

Set Shell "bash"
Set FontFamily "FiraCode Nerd Font"
Set FontSize 16
Set Width 1200
Set Height 480
Set LetterSpacing 0
Set LineHeight 1.2
Set Padding 20
Set Theme "Dracula"
Set MarginFill "#fe00fe"
Set WindowBar Colorful
Set Margin 20
Set BorderRadius 10

# Compress — note the before/after size, ratio, and throughput summary
Type "fqxv compress demo.fastq.gz -o demo.fqxv"
Enter
Sleep 3s

# Restore it — verified lossless on decode
Type "fqxv decompress demo.fqxv -o restored.fastq"
Enter
Sleep 3s
