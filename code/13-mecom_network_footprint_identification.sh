
# ----------------------------------------------------------------------------------------------------------------
# footprinting analysis
# ----------------------------------------------------------------------------------------------------------------
setwd("/Users/fyu/Documents/GitHub/")
# define function to generate footprint plot
plotFP <- function (Profile, Profile2="nodata", Mlen = 8, title="Footprint of FLI1 motif", xlab = "Dist. to motif (bp)", ylab = "Cut-site probability", 
    legTitle = "", legPos = "topright") 
{
    par(cex = 1)
    S <- length(Profile)
    W <- ((S/2) - Mlen)/2
    plot((1:(S/2)), apply(data.frame(Profile[1:(S/2)], Profile[(S/2 + 1):S]), 1, mean), axes = "F", xlab = "", 
        ylab = "", ylim = c(0, max(Profile) * 1.12), type = "l", 
        lwd = 2, col = "darkred", main=title)
    # lines((1:(S/2)), Profile[(S/2 + 1):S], lwd = 2, type = "l", 
    #     col = "darkred")
    if(Profile2!="nodata"){
      lines((1:(S/2)), apply(data.frame(Profile2[1:(S/2)], Profile2[(S/2 + 1):S]), 1, mean), ylim = c(0, max(Profile) * 1.12), 
        lwd = 2, col = "lightgray")
    }

    axis(1, at = seq(1, W, len = 3), labels = -(W + 1 - seq(1, 
        W + 1, len = 3)), padj = -1, tck = -0.01)
    axis(1, at = W + Mlen + seq(1, W, len = 3), labels = seq(0, 
        W, len = 3), padj = -1, tck = -0.01)
    axis(1, at = W + seq(1, Mlen), labels = NA, padj = -1.2, 
        tck = +0.01, col = "purple4")
    axis(2, padj = 1, tck = -0.02)
    abline(v = c(W, W + Mlen + 1), lty = 2)
    mtext(xlab, 1, padj = 1.6, cex = 1)
    mtext(ylab, 2, padj = -2, cex = 1)
    # if (missing(legTitle)) 
    #     legend(legPos, c("For. strand", "Rev. strand"), lwd = 2, 
    #         lty = c(1, 1), col = c("darkblue", "darkred"))
    # else {
    #     legend(legPos, legTitle, bty = "n")
    # }
}


### 
### down-regulated 
### 

# find the De novo motifs from meme analysis
# 1.(fimo_out_1) ETS
# 2.(fimo_out_2) CTCF
# 4.AACCACA (fimo_out_5) RUNX
# 5.TGACTCAG (fimo_out_6) JUN
# 6.CCCCGCCC (fimo_out_8) KLF

workdir=mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_down
TF_names="ETS CTCF RUNX JUN KLF"
data_dir=mecom_var_sankaran/data/motif_denovo/down_mast0035_interactions_10x_hocomoco_meme-chip
gff2bed < $data_dir/fimo_out_1/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $workdir/ETS_fimo.bed
gff2bed < $data_dir/fimo_out_2/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $workdir/CTCF_fimo.bed
gff2bed < $data_dir/fimo_out_5/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $workdir/RUNX_fimo.bed
gff2bed < $data_dir/fimo_out_6/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $workdir/JUN_fimo.bed
gff2bed < $data_dir/fimo_out_8/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $workdir/KLF_fimo.bed
hsc_bam=mecom_var_sankaran/data/footprint_analysis/HSC-1.merge.bam # Note that the bam file is too large to upload, download this file to the assigned folder first by yourself
hsc_dir=$workdir/hsc; mkdir -p $hsc_dir
centipede_r_script=mecom_var_sankaran/data/footprint_analysis/fimo_0005/run_centipede_parker.R
cd $hsc_dir
# ETS
make_cut_matrix -v -b '(36-149 150-324 325-400 1)' -d -o 0 -r 254 -p 1 -f 3 -F 4 -F 8 -q 0 $hsc_bam $workdir/ETS_fimo.bed > $workdir/hsc_ETS_fimo.cuts.freq.txt
Rscript $centipede_r_script $workdir/hsc_ETS_fimo.cuts.freq.txt $workdir/ETS_fimo.bed $hsc_dir/hsc_ETS_fimo.pdf 8
setwd("/Users/fyu/Library/Mobile Documents/com~apple~CloudDocs/Documents/Vijay/project/09022020-MECOM_Richard/data/02262021-downMast_re-analysis/fp/meme/mast0035_down")
ETS_cp_hsc <- read.table("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_down/hsc/hsc_ETS_fimo.pdf.lambda.txt")[, 1]
pdf("meme_500_mast0035_down_ETS_fp.pdf")
plotFP(ETS_cp_hsc, title="Footprint of ETS motif", Mlen = 8)
dev.off()
# CTCF
make_cut_matrix -v -b '(36-149 150-324 325-400 1)' -d -o 0 -r 258 -p 1 -f 3 -F 4 -F 8 -q 0 $hsc_bam $workdir/CTCF_fimo.bed > $workdir/hsc_CTCF_fimo.cuts.freq.txt
Rscript $centipede_r_script $workdir/hsc_CTCF_fimo.cuts.freq.txt $workdir/CTCF_fimo.bed $hsc_dir/hsc_CTCF_fimo.pdf 8
CTCF_cp_hsc <- read.table("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_down/hsc/hsc_CTCF_fimo.pdf.lambda.txt")[, 1]
pdf("meme_500_mast0035_down_CTCF_fp.pdf")
plotFP(CTCF_cp_hsc, title="Footprint of CTCF motif", Mlen = 16)
dev.off()
# RUNX
make_cut_matrix -v -b '(36-149 150-324 325-400 1)' -d -o 0 -r 254 -p 1 -f 3 -F 4 -F 8 -q 0 $hsc_bam $workdir/RUNX_fimo.bed > $workdir/hsc_RUNX_fimo.cuts.freq.txt
Rscript $centipede_r_script $workdir/hsc_RUNX_fimo.cuts.freq.txt $workdir/RUNX_fimo.bed $hsc_dir/hsc_RUNX_fimo.pdf 8
RUNX_cp_hsc <- read.table("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_down/hsc/hsc_RUNX_fimo.pdf.lambda.txt")[, 1]
pdf("meme_500_mast0035_down_RUNX_fp.pdf")
plotFP(RUNX_cp_hsc, title="Footprint of RUNX motif", Mlen = 8)
dev.off()
# JUN
make_cut_matrix -v -b '(36-149 150-324 325-400 1)' -d -o 0 -r 254 -p 1 -f 3 -F 4 -F 8 -q 0 $hsc_bam $workdir/JUN_fimo.bed > $workdir/hsc_JUN_fimo.cuts.freq.txt
Rscript $centipede_r_script $workdir/hsc_JUN_fimo.cuts.freq.txt $workdir/JUN_fimo.bed $hsc_dir/hsc_JUN_fimo.pdf 8
JUN_cp_hsc <- read.table("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_down/hsc/hsc_JUN_fimo.pdf.lambda.txt")[, 1]
pdf("meme_500_mast0035_down_JUN_fp.pdf")
plotFP(JUN_cp_hsc, title="Footprint of JUN motif", Mlen = 8)
dev.off()
# KLF
make_cut_matrix -v -b '(36-149 150-324 325-400 1)' -d -o 0 -r 254 -p 1 -f 3 -F 4 -F 8 -q 0 $hsc_bam $workdir/KLF_fimo.bed > $workdir/hsc_KLF_fimo.cuts.freq.txt
Rscript $centipede_r_script $workdir/hsc_KLF_fimo.cuts.freq.txt $workdir/KLF_fimo.bed $hsc_dir/hsc_KLF_fimo.pdf 8
KLF_cp_hsc <- read.table("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_down/hsc/hsc_KLF_fimo.pdf.lambda.txt")[, 1]
pdf("meme_500_mast0035_down_KLF_fp.pdf")
plotFP(KLF_cp_hsc, title="Footprint of KLF motif", Mlen = 8)
dev.off()


### 
### up-regulated 
### 

# find the De novo motifs from meme analysis
# 1.(fimo_out_1) ETS
# 2.(fimo_out_2) CTCF
# 3.(fimo_out_3) RUNX
# 5.TGACTCAG (fimo_out_6) JUN
# 6.CCCCGCCC (fimo_out_8) KLF

workdir=mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_up; mkdir -p $workdir
TF_names="ETS CTCF SP RUNX JUN KLF"
data_dir=mecom_var_sankaran/data/motif_denovo/up_mast0035_interactions_10x_hocomoco_meme-chip 
gff2bed < $data_dir/fimo_out_1/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $workdir/ETS_fimo.bed
gff2bed < $data_dir/fimo_out_2/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $workdir/CTCF_fimo.bed
gff2bed < $data_dir/fimo_out_3/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $workdir/RUNX_fimo.bed
gff2bed < $data_dir/fimo_out_6/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $workdir/JUN_fimo.bed
gff2bed < $data_dir/fimo_out_8/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $workdir/KLF_fimo.bed
# gff2bed < $data_dir/fimo_out_10/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $workdir/GATA_fimo.bed
hsc_bam=mecom_var_sankaran/data/footprint_analysis/HSC-1.merge.bam
hsc_dir=$workdir/hsc; mkdir -p $hsc_dir
centipede_r_script=mecom_var_sankaran/data/footprint_analysis/fimo_0005/run_centipede_parker.Rcd $hsc_dir
# ETS
make_cut_matrix -v -b '(36-149 150-324 325-400 1)' -d -o 0 -r 254 -p 1 -f 3 -F 4 -F 8 -q 0 $hsc_bam $workdir/ETS_fimo.bed > $workdir/hsc_ETS_fimo.cuts.freq.txt
Rscript $centipede_r_script $workdir/hsc_ETS_fimo.cuts.freq.txt $workdir/ETS_fimo.bed $hsc_dir/hsc_ETS_fimo.pdf 8
# R
setwd("/Users/fyu/Library/Mobile Documents/com~apple~CloudDocs/Documents/Vijay/project/09022020-MECOM_Richard/data/02262021-downMast_re-analysis/fp/meme/mast0035_up")
ETS_cp_hsc <- read.table("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_up/hsc/hsc_ETS_fimo.pdf.lambda.txt")[, 1]
pdf("meme_500_mast0035_up_ETS_fp.pdf")
plotFP(ETS_cp_hsc, title="Footprint of ETS motif", Mlen = 8)
dev.off()

# CTCF
make_cut_matrix -v -b '(36-149 150-324 325-400 1)' -d -o 0 -r 258 -p 1 -f 3 -F 4 -F 8 -q 0 $hsc_bam $workdir/CTCF_fimo.bed > $workdir/hsc_CTCF_fimo.cuts.freq.txt
Rscript $centipede_r_script $workdir/hsc_CTCF_fimo.cuts.freq.txt $workdir/CTCF_fimo.bed $hsc_dir/hsc_CTCF_fimo.pdf 8
# R
CTCF_cp_hsc <- read.table("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_up/hsc/hsc_CTCF_fimo.pdf.lambda.txt")[, 1]
pdf("meme_500_mast0035_up_CTCF_fp.pdf")
plotFP(CTCF_cp_hsc, title="Footprint of CTCF motif", Mlen = 16)
dev.off()

# RUNX
make_cut_matrix -v -b '(36-149 150-324 325-400 1)' -d -o 0 -r 254 -p 1 -f 3 -F 4 -F 8 -q 0 $hsc_bam $workdir/RUNX_fimo.bed > $workdir/hsc_RUNX_fimo.cuts.freq.txt
Rscript $centipede_r_script $workdir/hsc_RUNX_fimo.cuts.freq.txt $workdir/RUNX_fimo.bed $hsc_dir/hsc_RUNX_fimo.pdf 8
# R
RUNX_cp_hsc <- read.table("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_up/hsc/hsc_RUNX_fimo.pdf.lambda.txt")[, 1]
pdf("meme_500_mast0035_up_RUNX_fp.pdf")
plotFP(RUNX_cp_hsc, title="Footprint of RUNX motif", Mlen = 8)
dev.off()

# JUN
make_cut_matrix -v -b '(36-149 150-324 325-400 1)' -d -o 0 -r 254 -p 1 -f 3 -F 4 -F 8 -q 0 $hsc_bam $workdir/JUN_fimo.bed > $workdir/hsc_JUN_fimo.cuts.freq.txt
Rscript $centipede_r_script $workdir/hsc_JUN_fimo.cuts.freq.txt $workdir/JUN_fimo.bed $hsc_dir/hsc_JUN_fimo.pdf 8
# R
JUN_cp_hsc <- read.table("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_up/hsc/hsc_JUN_fimo.pdf.lambda.txt")[, 1]
pdf("meme_500_mast0035_up_JUN_fp.pdf")
plotFP(JUN_cp_hsc, title="Footprint of JUN motif", Mlen = 8)
dev.off()

# KLF
make_cut_matrix -v -b '(36-149 150-324 325-400 1)' -d -o 0 -r 254 -p 1 -f 3 -F 4 -F 8 -q 0 $hsc_bam $workdir/KLF_fimo.bed > $workdir/hsc_KLF_fimo.cuts.freq.txt
Rscript $centipede_r_script $workdir/hsc_KLF_fimo.cuts.freq.txt $workdir/KLF_fimo.bed $hsc_dir/hsc_KLF_fimo.pdf 8
# R
KLF_cp_hsc <- read.table("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_up/hsc/hsc_KLF_fimo.pdf.lambda.txt")[, 1]
pdf("meme_500_mast0035_up_KLF_fp.pdf")
plotFP(KLF_cp_hsc, title="Footprint of KLF motif", Mlen = 8)
dev.off()


### 
### all -deg
### 

# find the De novo motifs from meme analysis
# 1.(fimo_out_1) ETS
# 2.(fimo_out_6) CTCF
# 3.(fimo_out_2) RUNX
# 5.TGACTCAG (fimo_out_3) JUN
# 6.CCCCGCCC (fimo_out_4) KLF
# 7. (fimo_out_12) GATA

workdir=mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_all; mkdir -p $workdir
TF_names="ETS CTCF SP RUNX JUN KLF"
data_dir=/Users/fyu/Documents/GitHub/mecom_var_sankaran/data/motif_denovo/all_mast0035_interactions_10x_hocomoco_meme-chip 
gff2bed < $data_dir/fimo_out_1/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $workdir/ETS_fimo.bed
gff2bed < $data_dir/fimo_out_6/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $workdir/CTCF_fimo.bed
gff2bed < $data_dir/fimo_out_2/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $workdir/RUNX_fimo.bed
gff2bed < $data_dir/fimo_out_3/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $workdir/JUN_fimo.bed
gff2bed < $data_dir/fimo_out_4/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $workdir/KLF_fimo.bed
# /broad/software/free/Linux/redhat_7_x86_64/pkgs/meme_4.11.3_1/bin/fimo --parse-genomic-coord --thresh 5e-4 --verbosity 1 --oc /broad/sankaranlab/fyu/HemeMap/02262021-downMast_re-analysis/motif_denovo/all_mast0035_interactions_10x_hocomoco_meme-chip/fimo_out_12_0005 --bgfile /broad/sankaranlab/fyu/HemeMap/02262021-downMast_re-analysis/motif_denovo/all_mast0035_interactions_10x_hocomoco_meme-chip/background --motif AGAWAA /broad/sankaranlab/fyu/HemeMap/02262021-downMast_re-analysis/motif_denovo/all_mast0035_interactions_10x_hocomoco_meme-chip/dreme_out/dreme.xml /broad/sankaranlab/fyu/HemeMap/02262021-downMast_re-analysis/motif_denovo/all_mast0035_interactions_10x_hocomoco_meme-chip/20210227-deg_interactions_peaks_all_mast0035_fa_1based.fa
gff2bed < $data_dir/fimo_out_12_0005/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $workdir/GATA_fimo.bed
# gff2bed < $data_dir/fimo_out_10/fimo.gff | awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}' > $workdir/GATA_fimo.bed
hsc_bam=mecom_var_sankaran/data/footprint_analysis/HSC-1.merge.bam
hsc_dir=$workdir/hsc; mkdir -p $hsc_dir
centipede_r_script=mecom_var_sankaran/data/footprint_analysis/fimo_0005/run_centipede_parker.R
cd $hsc_dir
# ETS
make_cut_matrix -v -b '(36-149 150-324 325-400 1)' -d -o 0 -r 254 -p 1 -f 3 -F 4 -F 8 -q 0 $hsc_bam $workdir/ETS_fimo.bed > $workdir/hsc_ETS_fimo.cuts.freq.txt
Rscript $centipede_r_script $workdir/hsc_ETS_fimo.cuts.freq.txt $workdir/ETS_fimo.bed $hsc_dir/hsc_ETS_fimo.pdf 8
# R
setwd("/Users/fyu/Library/Mobile Documents/com~apple~CloudDocs/Documents/Vijay/project/09022020-MECOM_Richard/data/02262021-downMast_re-analysis/fp/meme/mast0035_all")
ETS_cp_hsc <- read.table("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_all/hsc/hsc_ETS_fimo.pdf.lambda.txt")[, 1]
pdf("meme_500_mast0035_all_ETS_fp.pdf")
plotFP(ETS_cp_hsc, title="Footprint of ETS motif", Mlen = 8)
dev.off()

# CTCF
make_cut_matrix -v -b '(36-149 150-324 325-400 1)' -d -o 0 -r 258 -p 1 -f 3 -F 4 -F 8 -q 0 $hsc_bam $workdir/CTCF_fimo.bed > $workdir/hsc_CTCF_fimo.cuts.freq.txt
Rscript $centipede_r_script $workdir/hsc_CTCF_fimo.cuts.freq.txt $workdir/CTCF_fimo.bed $hsc_dir/hsc_CTCF_fimo.pdf 8
# R
CTCF_cp_hsc <- read.table("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_all/hsc/hsc_CTCF_fimo.pdf.lambda.txt")[, 1]
pdf("meme_500_mast0035_all_CTCF_fp.pdf")
plotFP(CTCF_cp_hsc, title="Footprint of CTCF motif", Mlen = 16)
dev.off()

# RUNX
make_cut_matrix -v -b '(36-149 150-324 325-400 1)' -d -o 0 -r 254 -p 1 -f 3 -F 4 -F 8 -q 0 $hsc_bam $workdir/RUNX_fimo.bed > $workdir/hsc_RUNX_fimo.cuts.freq.txt
Rscript $centipede_r_script $workdir/hsc_RUNX_fimo.cuts.freq.txt $workdir/RUNX_fimo.bed $hsc_dir/hsc_RUNX_fimo.pdf 8
# R
RUNX_cp_hsc <- read.table("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_all/hsc/hsc_RUNX_fimo.pdf.lambda.txt")[, 1]
pdf("meme_500_mast0035_all_RUNX_fp.pdf")
plotFP(RUNX_cp_hsc, title="Footprint of RUNX motif", Mlen = 8)
dev.off()

# JUN
make_cut_matrix -v -b '(36-149 150-324 325-400 1)' -d -o 0 -r 254 -p 1 -f 3 -F 4 -F 8 -q 0 $hsc_bam $workdir/JUN_fimo.bed > $workdir/hsc_JUN_fimo.cuts.freq.txt
Rscript $centipede_r_script $workdir/hsc_JUN_fimo.cuts.freq.txt $workdir/JUN_fimo.bed $hsc_dir/hsc_JUN_fimo.pdf 8
# R
JUN_cp_hsc <- read.table("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_all/hsc/hsc_JUN_fimo.pdf.lambda.txt")[, 1]
pdf("meme_500_mast0035_all_JUN_fp.pdf")
plotFP(JUN_cp_hsc, title="Footprint of JUN motif", Mlen = 8)
dev.off()

# KLF
make_cut_matrix -v -b '(36-149 150-324 325-400 1)' -d -o 0 -r 254 -p 1 -f 3 -F 4 -F 8 -q 0 $hsc_bam $workdir/KLF_fimo.bed > $workdir/hsc_KLF_fimo.cuts.freq.txt
Rscript $centipede_r_script $workdir/hsc_KLF_fimo.cuts.freq.txt $workdir/KLF_fimo.bed $hsc_dir/hsc_KLF_fimo.pdf 8
# R
KLF_cp_hsc <- read.table("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_all/hsc/hsc_KLF_fimo.pdf.lambda.txt")[, 1]
pdf("meme_500_mast0035_all_KLF_fp.pdf")
plotFP(KLF_cp_hsc, title="Footprint of KLF motif", Mlen = 8)
dev.off()

# GATA
make_cut_matrix -v -b '(36-149 150-324 325-400 1)' -d -o 0 -r 254 -p 1 -f 3 -F 4 -F 8 -q 0 $hsc_bam $workdir/GATA_fimo.bed > $workdir/hsc_GATA_fimo.cuts.freq.txt
Rscript $centipede_r_script $workdir/hsc_GATA_fimo.cuts.freq.txt $workdir/GATA_fimo.bed $hsc_dir/hsc_GATA_fimo.pdf 8
# R
GATA_cp_hsc <- read.table("mecom_var_sankaran/data/footprint_analysis/fimo_0005/mast0035_all/hsc/hsc_GATA_fimo.pdf.lambda.txt")[, 1]
pdf("meme_500_mast0035_all_GATA_fp.pdf")
plotFP(GATA_cp_hsc, title="Footprint of GATA motif", Mlen = 8)
dev.off()

