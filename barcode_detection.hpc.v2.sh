# v2: this is the tidied up and later version of 
# /Users/lany/mount_drives/HPC_Projects/0/0ddcec61-75f4-4347-b696-e96ec587183c/data/results/YL/PIP-3583/Analysis/Barcode_analysis/run2/sbatch_barcode_detection.hpc.v2.sh

projects_dir=/ei/projects/0/0ddcec61-75f4-4347-b696-e96ec587183c
analysis_dir=$projects_dir/data/results/YL/PIP-3583/Analysis/Barcode_analysis
plateid=plate6
plateid_old=plate6fullLIBk562

# ===== 1. merge and subsample =====
if true; then
    # make symlink of the raw data
    if true; then
        rawdir=$projects_dir/data/raw/PIP-3583/250320_VL00510_81_AAGTGTFM5
        symlinkdir=$projects_dir/data/results/YL/PIP-3583/Reads/Raw
        cd $rawdir
        for fastq in `ls *_R1.fastq.gz`; do
            newname=`echo $fastq | cut -d"_" -f 3`
            ln -sfv $rawdir/$fastq $symlinkdir/$newname"_R1.fastq.gz"
        done 
        for fastq in `ls *_R2.fastq.gz`; do
            newname=`echo $fastq | cut -d"_" -f 3`
            ln -svf $rawdir/$fastq $symlinkdir/$newname"_R2.fastq.gz"
        done 
    fi 

    # merge R1 and R2
    if true; then
        merged_dir=$projects_dir/data/results/YL/PIP-3583/Reads/Raw_merged/
        cd $analysis_dir
        symlinkdir=$projects_dir/data/results/YL/PIP-3583/Reads/Raw

        for plateid in plate2 plate3 plate4 plate5 plate6fullLIBk562; do
        # source bbmap-39.13_CBG
            sbatch -p ei-short -J bbmerge_${plateid} -o bbmerge_${plateid}.%j.%N.log -c 1 --mem 20G \
                --mail-type=ALL --mail-user=yuxuan.lan@earlham.ac.uk \
                --wrap "source bbmap-39.13_CBG && cd $analysis_dir && \
                    bbmerge.sh in1=$symlinkdir/${plateid}_R1.fastq.gz \
                        in2=$symlinkdir/${plateid}_R2.fastq.gz \
                        out=$merged_dir/${plateid}_merged.fastq.gz \
                        outu1=$merged_dir/${plateid}_R1_unmerged.fastq.gz \
                        outu2=$merged_dir/${plateid}_R2_unmerged.fastq.gz" && sleep 30
            
        done
    fi 

fi
