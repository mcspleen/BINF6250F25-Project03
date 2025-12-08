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

# Personal Reflections

## Group Leader
This was my first encounter with Gibbs Sampling, and it is a concept I will probably continue struggling with at least for a couple more weeks. I get dizzy with the several loops within loops we had to do for this algorithm and keeping track of what it all means. Tiange definitely carried me when mapping out the direction we needed to take to complete the GibbsMotifFinder function, but when focusing in on writing out the individual lines of code I was thankfully able to contribute ideas and have an understanding of what was at least happening in that block of code. I am hoping that in a few weeks I can revisit this project with confidence and make any changes that we missed in this first rendition of our Gibbs Sampling algorithm.

## Other member - Tiange
This project was more challenging than the previous two. However, through some research and productive discussions with Jason, I gradually pinpointed the core solution: meticulous comprehension of the existing code modules (including their inputs, functions, and outputs) and establishing the logical flow between them. Before proceeding to implement the GibbsMotifFinder, I spent a considerable amount of time configuring the environment and importing necessary modules, some of which required sourcing from specific channels. This process provided valuable, hands-on exposure to environment management and configuration, which is an essential complement to the algorithmic, statistical, and IDE skills practiced during the assignment itself.

# Generative AI Appendix
1.  Googled "yield" (python)\
    [**`yield` pauses and saves state; `return` terminates.**](https://www.google.com/search?client=safari&sca_esv=861dc178219c5a10&rls=en&q=yield+pauses+and+saves+state%3B+return+terminates.&sa=X&ved=2ahUKEwiijL7t9PuPAxXDg4kEHahPOgQQxccNegUIiwEQAw&mstk=AUtExfBHwhHoWRb-HlKNjTlQ8FFThnClkWMlyCIYfbtR6UqIVMYwpyO_7Yh-letMs4CI9GgvoM5MWi85nvdbCu1zobFviPPOl2J37uByYeez1fO407msxk0BFSjJwPLK2Ti4Yu9k9dazT-Ypi_vUsZ7adWMvYWw-Qo2Y8o1oZw3BSMm33f4&csui=3)

    A `return` statement immediately exits the function, discarding its state. `yield` allows the function to be paused and resumed, enabling it to produce multiple values over time.

    [**`yield` creates a generator object; `return` returns a value.**](https://www.google.com/search?client=safari&sca_esv=861dc178219c5a10&rls=en&q=yield+creates+a+generator+object%3B+return+returns+a+value.&sa=X&ved=2ahUKEwiijL7t9PuPAxXDg4kEHahPOgQQxccNegUIjgEQAw&mstk=AUtExfBHwhHoWRb-HlKNjTlQ8FFThnClkWMlyCIYfbtR6UqIVMYwpyO_7Yh-letMs4CI9GgvoM5MWi85nvdbCu1zobFviPPOl2J37uByYeez1fO407msxk0BFSjJwPLK2Ti4Yu9k9dazT-Ypi_vUsZ7adWMvYWw-Qo2Y8o1oZw3BSMm33f4&csui=3)

    A function with `yield` becomes a generator function, which, when called, returns an iterable generator object. A function with `return` simply returns a single value.

2.  Gemini: save a list of selected sequences to a FASTA file

3.  Gemini: debug Test function\
    NameError: name 'COMPTRANS' is not defined. Did you mean: 'COMP_TRANS'?

4.  Gemini: random functions

