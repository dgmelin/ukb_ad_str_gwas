
# Example: 
# grep ADD replica.has_ad_icd10.glm.linear | awk -f format_to_meta.awk > tometa.txt
BEGIN{
    OFS="\t"; 
    FS="\t"; 
    print "CHR", "BP", "SNP", "A1", "A2", "OBS_CT", "BETA", "SE", "P"
}
$0~/^#/{next}
{
    # Target is CHR BP SNP A2 A1 OBS_CT BETA SE (Z_STAT) P -> (or similar)
    # mandatory is SNP, OR (or BETA), SE, P
    #skip non ADD lines
    if ($7!~/ADD/){next}
    # Source is #CHROM(1) POS(2) ID(3) REF(4) ALT(5) A1(6) TEST(7) OBS_CT(8) BETA(9) SE(10) T_STAT(11) P(12) ERRCODE(13)
    a1=$6
    a2=$4
    if (a2 == a1){
        a2=$5 
    }
    print $1, $2, $3, a1, a2, $8, $9, $10, $12
}