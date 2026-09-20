// test_parser.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "virovore.h"

// Tell the compiler that our evolutionary loop function exists in ga_loop.c
void run_evolutionary_loop(
    int generations,
    int pop_size,
    double mutation_rate,
    int tournament_size,
    unsigned int random_seed,
    int elite_count,
    double final_mutation_rate
);

int main(void) {
    printf("[*] Starting Ultra-Fidelity Silicon Virovore engine via test_parser.c...\n");
    
    // Explicit native-loop configuration; the mutation schedule is separate
    // from Python ARISE/ALE adaptation.
    run_evolutionary_loop(10000, 2000, 0.08, 7, 616, 200, 0.005);
    
    return 0;
}
