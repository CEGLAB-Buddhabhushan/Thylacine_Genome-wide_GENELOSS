for i in `cat L.transcripts.lst`; do ./Confirm_Loss.sh $i; mv $i L; done
