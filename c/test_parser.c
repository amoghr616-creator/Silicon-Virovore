// test_parser.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "virovore.h"

// Tell the compiler that our evolutionary loop function exists in ga_loop.c
void run_evolutionary_loop(int generations, int pop_size, double mutation_rate, int tournament_size);

int main(void) {
    printf("[*] Starting Ultra-Fidelity Silicon Virovore engine via test_parser.c...\n");
    
    // SCALE UP THE POPULATION: 
    // Run the genetic algorithm: 10,000 generations, 2,000 population size, 8% base mutation, tournament size 7
    run_evolutionary_loop(10000, 2000, 0.08, 7);
    
    return 0;
}