# Script to execute fimo for motif search in candidate mouse-specific enhancers and aligned non-functional human regions

# MEME version: meme-5.4.0

fimo --thresh 0.0001 --o ./fimo_out/hg38_non_functional_Cebpa  mouse_Cebpa_CISBP_2.00.meme hg38_non_functional_ov_mm10_CEBPA.fa

fimo --thresh 0.0001 --o ./fimo_out/mm10_enh_ov_CEBPA_Cebpa  mouse_Cebpa_CISBP_2.00.meme mm10_enh_ov_CEBPA.fa

fimo --thresh 0.0001 --o ./fimo_out/hg38_non_functional_Hnf4a  mouse_Hnf4a_CISBP_2.00.meme hg38_non_functional_ov_mm10_Hnf4a.fa

fimo --thresh 0.0001 --o ./fimo_out/mm10_enh_ov_Hnf4a_Hnf4a  mouse_Hnf4a_CISBP_2.00.meme mm10_enh_ov_Hnf4a.fa
