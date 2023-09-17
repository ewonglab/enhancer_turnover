# Script used for execute FIMO for motif search in candidate human-specific enhancers and the non-functional regions in the mouse genome that align to the candidate enhancers.

# meme version: meme-5.4.0

fimo --thresh 0.0001 --o ./fimo_out/mm10_non_functional_Cebpa  human_Cebpa_CISBP_2.00.meme mm10_non_functional_ov_hg38_CEBPA.fa

fimo --thresh 0.0001 --o ./fimo_out/hg38_enh_ov_CEBPA_Cebpa  human_Cebpa_CISBP_2.00.meme hg38_enh_ov_CEBPA.fa

fimo --thresh 0.0001 --o ./fimo_out/mm10_non_functional_Hnf4a  human_Hnf4a_CISBP_2.00.meme mm10_non_functional_ov_hg38_Hnf4a.fa

fimo --thresh 0.0001 --o ./fimo_out/hg38_enh_ov_Hnf4a_Hnf4a  human_Hnf4a_CISBP_2.00.meme hg38_enh_ov_Hnf4a.fa
