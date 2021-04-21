
# ----------------------------------------------------------------------------------------------------------------
# De Novo motif discovery using MEME suit
# ----------------------------------------------------------------------------------------------------------------
cd /Users/fyu/Documents/GitHub/
# ----------------------
# mecom down
workdir=mecom_var_sankaran/data/motif_denovo/downMast_re-analysis/motif_denovo ; mkdir -p $workdir
peak_data=mecom_var_sankaran/data/down0035_deg_info/down_mast0035_interactions_peak4.bed
genome_fa=/broad/sankaranlab/fyu/genome_data/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
peak_fa_name=deg_interactions_peaks_down_mast0035_fa
job_name=down_mast0035_interactions_10x
motif_data=/broad/sankaranlab/fyu/motif_data/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme
sort -k1,1V -k2,2n $peak_data | cut -f1-3 | bedtools getfasta -bed - -fi $genome_fa > $workdir/$peak_fa_name.fa
time awk -F ":" '{OFS=""; split($2,a,"-"); if(a[1]) print $1,":",a[1]+1,"-",a[2]; else print; }' $workdir/$peak_fa_name.fa > $workdir/${peak_fa_name}_1based.fa
meme-chip -oc $workdir/${job_name}_hocomoco_meme-chip -meme-p 4 -dreme-m 10 -meme-nmotifs 10 -db $motif_data $workdir/${peak_fa_name}_1based.fa

# mecom up
workdir=mecom_var_sankaran/data/motif_denovo/downMast_re-analysis/motif_denovo ; mkdir -p $workdir
peak_data=mecom_var_sankaran/data/down0035_deg_info/up_mast0035_interactions_peak4.bed
genome_fa=/broad/sankaranlab/fyu/genome_data/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
peak_fa_name=deg_interactions_peaks_up_mast0035_fa
job_name=up_mast0035_interactions_10x
motif_data=/broad/sankaranlab/fyu/motif_data/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme
sort -k1,1V -k2,2n $peak_data | cut -f1-3 | bedtools getfasta -bed - -fi $genome_fa > $workdir/$peak_fa_name.fa
time awk -F ":" '{OFS=""; split($2,a,"-"); if(a[1]) print $1,":",a[1]+1,"-",a[2]; else print; }' $workdir/$peak_fa_name.fa > $workdir/${peak_fa_name}_1based.fa
meme-chip -oc $workdir/${job_name}_hocomoco_meme-chip -meme-p 4 -dreme-m 10 -meme-nmotifs 10 -db $motif_data $workdir/${peak_fa_name}_1based.fa

# mecom all
workdir=mecom_var_sankaran/data/motif_denovo/downMast_re-analysis/motif_denovo ; mkdir -p $workdir
peak_data=mecom_var_sankaran/data/down0035_deg_info/all_mast0035_interactions_peak4.bed
genome_fa=/broad/sankaranlab/fyu/genome_data/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
peak_fa_name=deg_interactions_peaks_all_mast0035_fa
job_name=all_mast0035_interactions_10x
motif_data=/broad/sankaranlab/fyu/motif_data/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme
sort -k1,1V -k2,2n $peak_data | cut -f1-3 | bedtools getfasta -bed - -fi $genome_fa > $workdir/$peak_fa_name.fa
time awk -F ":" '{OFS=""; split($2,a,"-"); if(a[1]) print $1,":",a[1]+1,"-",a[2]; else print; }' $workdir/$peak_fa_name.fa > $workdir/${peak_fa_name}_1based.fa
meme-chip -oc $workdir/${job_name}_hocomoco_meme-chip -meme-p 4 -dreme-m 10 -meme-nmotifs 10 -db $motif_data $workdir/${peak_fa_name}_1based.fa

