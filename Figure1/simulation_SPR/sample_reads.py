import sys
import random
import itertools
import HTSeq
from itertools import izip

fraction = float( sys.argv[1] )
in1 = HTSeq.FastqReader( sys.argv[2] ) 
in2 = HTSeq.FastqReader( sys.argv[3] ) 
out1 = open( sys.argv[4], "w" )
out2 = open( sys.argv[5], "w" )
records=sum(1 for i in iter(in1))
rand_records =sorted(random.sample(range(0, records), int(fraction*records)))
for n, (read1, read2) in enumerate(izip(in1,in2)):
	if n in rand_records:
		read1.write_to_fastq_file( out1 )
		read1.write_to_fastq_file( out2 )
out1.close()
out2.close()

