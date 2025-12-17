BEGIN{
    OFS="\t";
    FS="\t";
    print "SNP", "A1", "A2", "freq", "b", "se", "p", "N"
}
# skip comments
$0~/^#/{next;}
NR==FNR{
    freq[$1]=$4; 
    next
}
{   
    if ($7 != "ADD"){
        next;
    }
    
    # Format is 
    # #CHROM	POS	ID	REF	ALT	A1	TEST	OBS_CT	BETA	SE	T_STAT	P	ERRCODE
    # Format should be
    # SNP A1 A2 freq b se p N

    # get frequencies
    split(freq[$3],fs,",")
    #take a1 from plink output
    a1=$6
    #check if it is ref or alt, assign a2 to other
    if (a1 == $4){
        # a1 == 'REF'
        a2=$5 
        a1_freq=fs[1]
    }else{
        # a1 == alt
        a2=$4
        a1_freq=fs[2]
    }
    print $3, a1, a2, a1_freq, $9, $10, $12, $8
}