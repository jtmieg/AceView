/^@/{next;}
{  # parse a SAM file generated by magicblast and list the introns and their support
    q = $3 ;     # q is the quey sequence, possibly a NC_number  chromosome
    pos = $4 ;   # first base of the alignment
    cigar = $6 ;
    k = 0 ;      # a place holder to parse integral numbers in decimal notation 
    da = -1 ;    # offset on the genome relative to pos, initialized left of first base
    # print "# " cigar;     # debugging trace
    for (i = 0 ; i < length (cigar) ; i++)
    {
	c = substr (cigar, i, 1) ;     # a piece of tobacco
	if (c >= "0" && c <= "9")      # is c a digit
	{
	    k = 10 * k + (c - 0) ;   # integrate it in k
	    continue ;
	}
	# print  "# ", k, c ;          # debugging trace
	type = c ;                     # c is a type
	dda = k ; k = 0 ;              # length of this feature
	if (type == "S")               # clipped, does not count on the genome
	    continue ;
 	if (type == "I")               # insertion, does not count on the genome
	    continue ;
        if (type == "N")               # register an intron
	{
	    a1 = pos + da + 1 ;        # first base of intron
	    a2 = pos + da + dda ;      # last base of intron

	    z = q "\t" a1 "\t" a2 ;    # u=intron full name
	    support[z] ++ ;            # integrate it in a hash table, awk style
	}
	da += dda                      # incorporate the feature length in the running offset
    }
}
END {                                  # export the cumultated supports
    print "# Target\tFirst base of intron\tLast base (i.e. the [g...g] of gt_ag)\tsupport"
    for (z in support)                       # loop on all introns
	printf ("%s\t%d\n", z, support[z]) ; # export the support
}
