#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "virovore.h"

// Forward declaration from engine.c
extern void evaluate_fitness(Variant *v, const char *target, const char *decoy1, const char *decoy2, void *params);

#define HERV_K_ENV "MKLAVDALLVTFAGSSDKKRR"
#define ALBUMIN "DAHKSEVAHRFKDLGEENFKALVL"
#define NCAM1 "MLQTKDLIWTLFFLGTAVS"

#define SEQ_LEN 21
static const char AA_ALPHABET[] = "ACDEFGHIKLMNPQRSTVWY";

static int is_valid_amino_acid(char c) {
    return (strchr(AA_ALPHABET, toupper((unsigned char)c)) != NULL);
}

/**
 * Executes full biophysical evaluation from engine.c on the given candidate.
 */
double c_check_sequence_fitness(const char* sequence) {
    if (sequence == NULL) return -999.0;
    
    int len = strlen(sequence);
    if (len < 9) return -1.0;
    
    for (int i = 0; i < len; i++) {
        if (!is_valid_amino_acid(sequence[i])) {
            return -1.0;
        }
    }
    
    // Instantiate temporary Variant structure
    Variant v = init_variant(len + 1);
    strcpy(v.sequence, sequence);
    v.length = len;

    // Call full biophysical evaluator in engine.c
    evaluate_fitness(&v, HERV_K_ENV, ALBUMIN, NCAM1, NULL);

    double score = v.fitness_score;

    // Clean up temporary variant buffer memory
    free(v.sequence);

    return score;
}

/**
 * C dynamic population generator
 */
void c_generate_mutated_population(const char* seed_sequence, char output_population[][SEQ_LEN + 1], int pop_size, double mutation_rate) {
    if (seed_sequence == NULL || strlen(seed_sequence) != SEQ_LEN) return;

    static int rand_seeded = 0;
    if (!rand_seeded) {
        srand((unsigned int)time(NULL));
        rand_seeded = 1;
    }

    for (int p = 0; p < pop_size; p++) {
        for (int i = 0; i < SEQ_LEN; i++) {
            double r = (double)rand() / RAND_MAX;
            if (r < mutation_rate) {
                output_population[p][i] = AA_ALPHABET[rand() % 20];
            } else {
                output_population[p][i] = seed_sequence[i];
            }
        }
        output_population[p][SEQ_LEN] = '\0';
    }
}

void c_generate_adaptive_population(
    const char *seed,
    double importance[],
    char output[][SEQ_LEN + 1],
    int pop_size)
{
    if (seed == NULL || strlen(seed) != SEQ_LEN) return;

    static int rand_seeded = 0;
    if (!rand_seeded) {
        srand((unsigned int)time(NULL));
        rand_seeded = 1;
    }

    for (int p = 0; p < pop_size; p++) {
        strcpy(output[p], seed);

        for (int i = 0; i < SEQ_LEN; i++) {
            /*
             * High importance
             * -> mutate less
             *
             * Low importance
             * -> mutate more
             */

            double mutation_rate = 0.05 + (1.0 - importance[i]) * 0.20;
            double r = (double)rand() / RAND_MAX;

            if (r < mutation_rate) {
                output[p][i] = AA_ALPHABET[rand() % 20];
            }
        }

        output[p][SEQ_LEN] = '\0';
    }
}