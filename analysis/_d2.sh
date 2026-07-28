#!/bin/bash
B="$1"; SAM="$2"
S=$(basename "$B" | sed 's/_combined.sorted.bam$//')
f(){ $SAM depth -a $1 -r "$2" "$B" 2>/dev/null | awk '{s+=$3;n++}END{if(n)printf "%.3f",s/n; else printf "NA"}'; }
a=$(f "-Q 20" chrY:1-4000000); b=$(f "" chrY:1-4000000)
c=$(f "-Q 20" chrY);           d=$(f "" chrY)
e=$(f "-Q 20" chr1:30000000-31000000)
printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$S" "$a" "$b" "$c" "$d" "$e"
