# Detect mid-file record dropout in a FASTQ: every intact record's "+" line
# (line 3 of 4) either repeats the header's read ID or is a bare "+". A header
# that names a different read than its "+" line proves records were dropped
# upstream, misaligning everything after the gap. Exit 1 if any mismatch found.
NR % 4 == 1 { h = substr($1, 2) }               # header line: strip leading @
NR % 4 == 3 {
    p = substr($1, 2)                            # "+" line: strip leading +
    if (p != "" && p != h) {
        print "MISALIGNED records in " f ": header @" h " but \"+\" line names +" p " (record " int(NR/4)+1 ")" > "/dev/stderr"
        bad = 1
    }
}
END { exit bad }