#!/bin/sh
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=80:00:00
#SBATCH --job-name=tree_analysis
#SBATCH -o launch_scripts/cluster_logs/atac.log 
#SBATCH --mem-per-cpu=5000MB

set -eo pipefail
module load R/4.0.5
TMPDIR=$(mktemp -d /tmp/$USER__$SLURM_JOBID-XXXX)

#Start the job here.
sbatch --wait launch_scripts/04.1-read_counts.sh
Rscript analysis/ATAC/MEGABULKS/0-extract_read_counts.R

Rscript analysis/ATAC/MEGABULKS/1.1-edgeR_megabulks.R
Rscript analysis/ATAC/MEGABULKS/1.2-edgeR_megabulks_plus_cn.R

Rscript analysis/ATAC/MEGABULKS/2.1-edgeR_megabulks_summary.R
Rscript analysis/ATAC/MEGABULKS/2.2-test_rnaseq_correlation.R
Rscript analysis/ATAC/MEGABULKS/2.3-edgeR_megabulks_plot_results.R

sbatch --wait launch_scripts/4-power_analysis_gather_data.sh

# SC data
Rscript analysis/ATAC/CLAIRE/Epigenetic/01-ATAC_regional_get_recurrent_peaks_to_test.R 
sbatch --wait analysis/ATAC/CLAIRE/Epigenetic/02-ATAC_regional_get_recurrent_peaks_readcounts.sh &
sbatch --wait analysis/ATAC/CLAIRE/Epigenetic/02-ATAC_regional_get_all_peaks_readcounts.sh &
wait
sbatch --wait analysis/ATAC/CLAIRE/Epigenetic/03-ATAC_regional_get_normfacs.sh
sbatch --wait analysis/ATAC/CLAIRE/Epigenetic/04-effect_of_purity_only_nfreads_LRT_runscript.sh &
sbatch --wait analysis/ATAC/CLAIRE/Epigenetic/04-effect_of_region_nfreads_LRT_runscript.sh &
sbatch --wait analysis/ATAC/CLAIRE/Epigenetic/04-effect_of_region_only_nfreads_LRT_runscript.sh &
sbatch --wait analysis/ATAC/CLAIRE/Epigenetic/04-effect_of_subtissues_for_adenoma.sh &
wait
sbatch --wait analysis/ATAC/CLAIRE/Epigenetic/05-ATAC_regional_get_recurrent_peak_coverage.sh 
Rscript analysis/ATAC/CLAIRE/Epigenetic/06-load_peak_track.R
Rscript analysis/ATAC/CLAIRE/Epigenetic/07-extract_deseq_results.R
Rscript analysis/ATAC/CLAIRE/Epigenetic/08-plot_results.R

echo "Tmp size is (after job run): $(du -sh "$TMPDIR")"
rm -rf "$TMPDIR"
