# Hidden Motif Discovery in Genomic Dataset

In this task, we create a genomic dataset with random sequence length and composition that contain a hidden motif with a certain mutation rate.

For the motif discovery, we check for the highest similarity motifs in a representative number of pairwise sequences (first few tens).
We assume that the motif has unknown length but ranges between 8-13 base pairs. Starting from minimum length=8 we iteratively search for all pairs of motifs of such length (kmers with k=8) between all the representative pairwise sequences. If any of those pairs of motifs has a similarity above a given threshold (say 0.8) both motifs are considered of high similarity.

Now, to establish the hidden motif, we will search across all the genomic dataset looking for those high similarity motifs. We discard the ones that are not present in all sequences with higher similarity than a given threshold. To choose the most probable motif of the ones remaining, we calculate the average similarity across all the sequences and we choose as a hidden motif the one with the highest average similarity among all. 

If more than one motif share the same average similarity or even slightly different, it is suggested to pick the highest length motif. Our motif generator do not nessacarily produce a hidden motif with the highest similarity . It is likely that a smaller sequence and usually a subsequence of our actual hidden motif can present a higher average similarity.
