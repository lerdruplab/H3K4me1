
# exctracting non-cg 
for sample in *.CX_report.txt
do
    awk '
    $6 != "CG" {
        cov = $4 + $5
        if (cov >= 5 && cov <= 500) {
            pct = 100 * $4 / cov
            print $0 "\t" pct
        }
    }' "$sample" > "${sample%.CX_report.txt}.nonCG.filtered.txt"
done


# calculating coverage in 10kb bins
for sample in *.nonCG.filtered.txt
do
    base=${sample%.nonCG.filtered.txt}

    awk 'BEGIN{OFS="\t"}
    {
        print $1,$2-1,$2,$4,$5
    }' "$sample" > "${base}.nonCG.counts.bed"

    bedtools map \
        -a /projects/lerdrup_group/data/k4me1/kubodata/binfiles/mm10_10KB_binning_sorted.bed \
        -b "${base}.nonCG.counts.bed" \
        -c 4,5 \
        -o sum,sum \
        > "${base}.nonCG.10kb.counts.bed"

    awk 'BEGIN{OFS="\t"}
    {
        meth=$4
        unmeth=$5

        if(meth=="." || unmeth=="." || meth+unmeth==0)
            pct="NA"
        else
            pct=100*meth/(meth+unmeth)

        print $1,$2,$3,meth,unmeth,pct
    }' "${base}.nonCG.10kb.counts.bed" \
    > "${base}.nonCG.filtered.10kbp.txt"

done
