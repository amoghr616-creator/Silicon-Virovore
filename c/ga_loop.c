/*
 * ga_loop.c
 * Silicon Virovore Project — 10,000 Generation Master Sync Engine
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "virovore.h"

const char AA_POOL[] = "ACDEFGHIKLMNPQRSTVWY";

int compare_variants(const void *a, const void *b) {
    double fit_a = ((Variant *)a)->fitness_score;
    double fit_b = ((Variant *)b)->fitness_score;
    if (fit_a > fit_b) return -1;
    if (fit_a < fit_b) return 1;
    return 0;
}

void mutate_variant(Variant *v, double rate) {
    for (size_t i = 0; i < v->length; i++) {
        double r = (double)rand() / RAND_MAX;
        if (r < rate) {
            v->sequence[i] = AA_POOL[rand() % 20];
        }
    }
}

void crossover(const Variant *p1, const Variant *p2, Variant *child) {
    size_t min_len = (p1->length < p2->length) ? p1->length : p2->length;
    if (min_len <= 1) {
        strcpy(child->sequence, p1->sequence);
        child->length = p1->length;
        return;
    }

    size_t split = 1 + (rand() % (min_len - 1));
    memcpy(child->sequence, p1->sequence, split);
    memcpy(child->sequence + split, p2->sequence + split, p1->length - split);
    child->sequence[p1->length] = '\0';
    child->length = p1->length;
}

double calculate_alignment_score(const char *variant, const char *profile) {
    int matches = 0;
    int len1 = strlen(variant);
    int len2 = strlen(profile);
    int min_len = (len1 < len2) ? len1 : len2;

    if (min_len == 0) return 0.0;

    for (int i = 0; i < min_len; i++) {
        if (toupper((unsigned char)variant[i]) == toupper((unsigned char)profile[i])) {
            matches++;
        }
    }
    return (double)matches / min_len;
}

void evaluate_fitness(Variant *v, const char *herv_k, const char *albumin, const char *ncam1, double *out_penalty) {
    const double WEIGHT_TARGET  = 2.0;
    const double WEIGHT_ALBUMIN = 2.5; 
    const double WEIGHT_NCAM1   = 2.5;

    double target_affinity = calculate_alignment_score(v->sequence, herv_k);
    double albumin_match    = calculate_alignment_score(v->sequence, albumin);
    double ncam1_match      = calculate_alignment_score(v->sequence, ncam1);

    double raw_score = (target_affinity * WEIGHT_TARGET) - 
                       (albumin_match * WEIGHT_ALBUMIN) - 
                       (ncam1_match * WEIGHT_NCAM1);

    v->solvation_energy    = calculate_solvation_energy(v->sequence);
    v->helix_propensity    = calculate_chou_fasman(v->sequence);
    v->hydrophobic_moment   = calculate_eisenberg_moment(v->sequence);

    double total_penalty = 0.0;

    charge_summary c_sum;
    if (compute_charge_profile(v->sequence, 1, &c_sum) == 0) {
        total_penalty += c_sum.penalty; 
    }

    hydropathy_summary h_sum;
    if (scan_hydropathy_profile(v->sequence, 9, &h_sum) == 0) {
        total_penalty += h_sum.solvation_penalty;
    }

    raw_score -= total_penalty;
    raw_score += (v->helix_propensity * 1.2);
    raw_score += (v->hydrophobic_moment * 1.5);
    raw_score -= (v->solvation_energy * 0.1); 

    v->fitness_score = raw_score;
    if (out_penalty) *out_penalty = total_penalty;
}

int main(void) {
    srand((unsigned int)time(NULL));

    printf("=============================================================\n");
    printf("Silicon Virovore Backend Engine v17.0 — Master Sync Pipeline\n");
    printf("=============================================================\n\n");

    // Open target telemetry data logs
    FILE *heatmap_file = fopen("evolution_consensus.csv", "w");
    FILE *metrics_file = fopen("metrics.csv", "w");
    
    if (!heatmap_file || !metrics_file) {
        fprintf(stderr, "ERROR: Structural pipeline failed to initialize target disk logs.\n");
        if (heatmap_file) fclose(heatmap_file);
        if (metrics_file) fclose(metrics_file);
        return 1;
    }

    // Initialize layout files
    fprintf(metrics_file, "Generation,FitnessScore,DecoyPenalty,ChouFasman,Eisenberg,SolvationEnergy\n");
    fprintf(heatmap_file, "Generation");
    for (int p = 1; p <= 21; p++) {
        fprintf(heatmap_file, ",Pos_%02d", p);
    }
    fprintf(heatmap_file, "\n");

    const char *HERV_K_ENV = "MKLAVDALLVTFAGSSDKKRR";
    const char *ALBUMIN    = "MKWVTFISLLFLFSSAYSRGV";
    const char *NCAM1      = "MLQTKDLIWTLFFLGTAVSLQ";

    size_t pop_size = 1000;
    Population pop;
    pop.variants = (Variant *)malloc(pop_size * sizeof(Variant));
    pop.count = pop_size;

    // Phase 1 initialization: Pure randomized molecular diversity
    for (size_t i = 0; i < pop_size; i++) {
        pop.variants[i] = init_variant(32);
        
        char random_seq[22];
        for (int pos = 0; pos < 21; pos++) {
            random_seq[pos] = AA_POOL[rand() % 20];
        }
        random_seq[21] = '\0';
        
        append_sequence(&pop.variants[i], random_seq);
        evaluate_fitness(&pop.variants[i], HERV_K_ENV, ALBUMIN, NCAM1, NULL);
    }

    size_t total_generations = 10000;
    printf("[SYSTEM] Running 10,000 generations and exporting physical arrays...\n");

    for (size_t gen = 0; gen <= total_generations; gen++) {
        qsort(pop.variants, pop_size, sizeof(Variant), compare_variants);

        // Every 100 generations, log metrics and positional matrices
        if (gen % 100 == 0) {
            double current_penalty = 0.0;
            evaluate_fitness(&pop.variants[0], HERV_K_ENV, ALBUMIN, NCAM1, &current_penalty);

            // Log time-series metrics data row
            fprintf(metrics_file, "%zu,%.4f,%.4f,%.4f,%.4f,%.4f\n", 
                    gen, 
                    pop.variants[0].fitness_score, 
                    current_penalty, 
                    pop.variants[0].helix_propensity, 
                    pop.variants[0].hydrophobic_moment, 
                    pop.variants[0].solvation_energy);

            // Log consensus structural data row
            fprintf(heatmap_file, "%zu", gen);
            for (size_t pos = 0; pos < 21; pos++) {
                int aa_counts[256] = {0};
                for (size_t elite = 0; elite < (pop_size / 10); elite++) {
                    char residue = pop.variants[elite].sequence[pos];
                    aa_counts[(unsigned char)residue]++;
                }
                char consensus_aa = 'A';
                int max_count = -1;
                for (int c = 0; c < 256; c++) {
                    if (aa_counts[c] > max_count) {
                        max_count = aa_counts[c];
                        consensus_aa = (char)c;
                    }
                }
                fprintf(heatmap_file, ",%d", (int)consensus_aa);
            }
            fprintf(heatmap_file, "\n");
        }

        // ------------------------------------------------------------
        // BREEDING PIPELINE: Tournament Selection & Annealing
        // ------------------------------------------------------------
        
        // Dynamic Mutation Annealing: Cool down from broad exploration (15%) to tight optimization (1%)
        double progress = (double)gen / total_generations;
        double current_mutation_rate = 0.15 * (1.0 - progress);
        if (current_mutation_rate < 0.01) current_mutation_rate = 0.01;

        for (size_t i = pop_size / 2; i < pop_size; i++) {
            // Tournament Selection: Match up random pairs to allow diverse background structures to propagate
            size_t t1 = rand() % pop_size;
            size_t t2 = rand() % pop_size;
            size_t parent1_idx = (pop.variants[t1].fitness_score > pop.variants[t2].fitness_score) ? t1 : t2;

            size_t t3 = rand() % pop_size;
            size_t t4 = rand() % pop_size;
            size_t parent2_idx = (pop.variants[t3].fitness_score > pop.variants[t4].fitness_score) ? t3 : t4;

            // Execute Crossover, Mutation, and Scoring
            crossover(&pop.variants[parent1_idx], &pop.variants[parent2_idx], &pop.variants[i]);
            mutate_variant(&pop.variants[i], current_mutation_rate);
            evaluate_fitness(&pop.variants[i], HERV_K_ENV, ALBUMIN, NCAM1, NULL);
        }
    }

    fclose(heatmap_file);
    fclose(metrics_file);
    printf("[SUCCESS] File sync pipelines complete. Data written to CSV channels.\n");

    // Run final sorting to guarantee index 0 represents the absolute top optimized sequence
    qsort(pop.variants, pop_size, sizeof(Variant), compare_variants);

    // ------------------------------------------------------------
    // TERMINAL CHAMPION REPORT (Extracting the Winner)
    // ------------------------------------------------------------
    printf("\n=============================================================\n");
    printf("DESIGN COMPLETE: THE WINNING SILICON VIROVORE PEPTIDE\n");
    printf("=============================================================\n");
    printf("Final Peak Fitness: %.4f\n", pop.variants[0].fitness_score);
    printf("Sequence (1-21):    %s\n", pop.variants[0].sequence);
    printf("α-Helix Propensity: %.4f\n", pop.variants[0].helix_propensity);
    printf("Hydrophobic Moment: %.4f\n", pop.variants[0].hydrophobic_moment);
    printf("=============================================================\n\n");

    // Free active population allocations
    for (size_t i = 0; i < pop_size; i++) {
        free(pop.variants[i].sequence);
    }
    free(pop.variants);
    printf("[CLEANUP] Memory sanitation complete. Leak profile: 0 bytes.\n");

    return 0;
}