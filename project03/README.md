# Introduction
This project implements the Gibbs Sampler algorithm in Python to discover a statistically significant sequence motif from a pre-filtered set of sequences. The sequences used are 50 base pair long promoter regions containing the Shine-Dalgarno sequence, extracted from the organism's FASTA and GFF files. The primary goal is to iteratively refine the motif's Position Frequency Matrix (PFM) until convergence, demonstrating an understanding of unsupervised motif finding.

# Pseudocode
```
FUNCTION GibbsMotifFinder(sequences, k, seed):
    SET N = number of sequences
    SET stability_window = 50
    SET convergence_threshold = 1e-5
    
    INITIALIZE motifs_init list by RANDOMLY selecting a k-mer from each sequence
    INITIALIZE ic_hist (Information Content history) list
    
    FOR samp_times FROM 0 TO 99:
        RANDOMLY select a sequence entry to 'leave out' (i)
        
        # 1. Build PFM/PWM from all motifs EXCEPT the one from sequence (i)
        SET current_motif = motifs_init without motif from index (i)
        SET pfm = build_pfm(current_motif, k)
        SET ic = pfm_ic(pfm)
        ic_hist.APPEND(ic)
        SET pwm = build_pwm(pfm)
        
        # 2. Resample the best motif for sequence (i)
        SET ith_entry = sequences[i]
        INITIALIZE kmer_score_list
        
        FOR kmer_index FROM 0 TO length(ith_entry) - k:
            SET motif_str = kmer at current index
            SET motif_score = score_kmer(motif_str, pwm)
            kmer_score_list.APPEND(motif_score)
            
        # Select the k-mer with the maximum score (best match to the new PWM)
        SET max_score = MAX(kmer_score_list)
        SET max_index = INDEX of max_score
        
        # 3. Update the motif list with the new, best-scoring k-mer
        motifs_init[i] = ith_entry[max_index : max_index + k]
        
        # 4. Check for Convergence (after 50 iterations)
        IF samp_times > 50:
            SET window = last stability_window IC values
            SET max_change = MAX(window) - MIN(window)
            
            IF max_change < convergence_threshold:
                PRINT "Convergence reached"
                RETURN pfm
                
        IF samp_times == 99:
            RETURN pfm
```

# Successes
Successfully implemented the core Gibbs Sampler logic: We successfully wrote the iterative core of the sampler, including the 'leave-one-out' PFM construction, k-mer scoring using the PWM, and updating the motif set based on the maximum score.

Integration of Existing Modules: We integrated pre-existing functions (build_pfm, build_pwm, score_kmer, pfm_ic) into the GibbsMotifFinder function, speeding up the implementation process.

Initial Convergence Check Implementation: We implemented a basic convergence check using Information Content (IC) stability over a 50-iteration window, ensuring the algorithm terminates when the motif stops significantly changing.

# Struggles
External Package Management: We experienced difficulties installing the seqlogo and pdf2svg packages correctly using reticulate::conda_install, particularly related to specifying the correct channels and resolving dependency conflicts with the forge=TRUE argument and Conda solver (requiring the addition of --solver=libmamba). Our program needs ipython environment to generate the result plot, where we had problems to locate and debug.



# Generative AI Appendix
The following generative AI resources were used for debugging, understanding concepts, and improving code efficiency:

Google Search/Snippets: Clarification on the fundamental difference between the Python commands yield (pauses and saves state; creates a generator object) and return (immediately exits; returns a single value).

Gemini (LLM): Assisted with debugging a NameError related to case sensitivity (COMPTRANS vs. COMP_TRANS), provided basic Python random function usage, and offered guidance on saving selected sequences to a FASTA file.
