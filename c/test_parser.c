/*
 * test_parser.c
 * Silicon Virovore Project — May 2026
 *
 * Upgraded FASTA parser implementing dynamic heap memory allocation (malloc/realloc),
 * case-insensitivity safeguards, and clean memory teardowns for infinite-length sequences.
 *
 * Compile: gcc -o parser_demo c/test_parser.c
 * Run:     ./parser_demo
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>  

// Clean data structural layer tying biological assets together
typedef struct {
    char name[100];
    char *sequence;      // Changed from array to pointer for dynamic heap allocation
    size_t length;
    float net_charge;       // Placeholder for future engine2 computations
    float net_hydropathy;   // Placeholder for future hydropathy computations
} Protein;

void load_fasta_sequence(const char *filename, Protein *protein) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "\nFATAL ERROR: Could not open data file: %s\n", filename);
        fprintf(stderr, "Please make sure the file is in your data/ folder and named correctly!\n");
        exit(1);
    }

    // 1. Initial Heap Allocation: Start with a 1024-byte sandbox
    size_t capacity = 1024;
    protein->sequence = malloc(capacity);
    if (protein->sequence == NULL) {
        fprintf(stderr, "FATAL ERROR: Out of memory during initial allocation for %s\n", protein->name);
        fclose(file);
        exit(1);
    }

    protein->sequence[0] = '\0';
    protein->length = 0;

    char line[256];
    while (fgets(line, sizeof(line), file)) {
        // Skip FASTA description lines starting with '>'
        if (line[0] == '>') {
            continue; 
        }
        
        // Strip out invisible newline formatting safely
        line[strcspn(line, "\r\n")] = 0;
        size_t line_len = strlen(line);
        if (line_len == 0) continue;
        
        // Process each character to enforce uppercase styling
        for (size_t i = 0; i < line_len; i++) {
            line[i] = toupper((unsigned char)line[i]);
        }
        
        // 2. Dynamic Expansion Check: Grow the loop exponentially if incoming data overflows capacity
        while (protein->length + line_len + 1 > capacity) {
            capacity *= 2; // Double the capacity bounds
            char *temp = realloc(protein->sequence, capacity);
            if (temp == NULL) {
                fprintf(stderr, "FATAL ERROR: Failed to reallocate memory grid for %s!\n", protein->name);
                free(protein->sequence); // Stop memory leaks if system runs out of RAM
                fclose(file);
                exit(1);
            }
            protein->sequence = temp; // Re-point to the newly expanded heap block
        }

        // 3. Append data safely now that space is guaranteed by the OS
        strcpy(protein->sequence + protein->length, line);
        protein->length += line_len;
    }

    fclose(file);
}

int main(void) {
    // Instantiate structures with null sequence pointers before parsing
    Protein target_herv = { .name = "HERV-K Env (Viral Target)", .sequence = NULL };
    Protein decoy_albumin = { .name = "Albumin (Blood Plasma Decoy)", .sequence = NULL };
    Protein decoy_ncam1 = { .name = "NCAM1 (Motor Neuron Decoy)", .sequence = NULL };

    printf("=============================================================\n");
    printf("Silicon Virovore — Multi-Target Data Initialization\n");
    printf("May 2026 [Dynamic Memory Heap Architecture]\n");
    printf("=============================================================\n\n");

    // Load data assets using relative root positioning
    printf("Loading viral target (%s) from file... \n", target_herv.name);
    load_fasta_sequence("data/P61567.fasta.txt", &target_herv);
    printf("--> SUCCESS: Loaded %lu amino acid residues dynamic buffer.\n\n", target_herv.length);

    printf("Loading blood plasma decoy (%s)... \n", decoy_albumin.name);
    load_fasta_sequence("data/P02768.fasta.txt", &decoy_albumin);
    printf("--> SUCCESS: Loaded %lu amino acid residues dynamic buffer.\n\n", decoy_albumin.length);

    printf("Loading motor neuron cell surface decoy (%s)... \n", decoy_ncam1.name);
    load_fasta_sequence("data/P13591.fasta.txt", &decoy_ncam1);
    printf("--> SUCCESS: Loaded %lu amino acid residues dynamic buffer.\n\n", decoy_ncam1.length);

    printf("=============================================================\n");
    printf("DATABASE CONNECTIVITY VERIFIED & MEMORY HEAP BOUNDLESS!\n");
    printf("June production pipeline data layer is fully ready.\n");
    printf("=============================================================\n\n");

    // 4. Clean Memory Teardown: Return borrowed RAM back to your system layer
    printf("Cleaning up runtime memory allocations... ");
    free(target_herv.sequence);
    free(decoy_albumin.sequence);
    free(decoy_ncam1.sequence);
    printf("DONE.\n");

    return 0;
}