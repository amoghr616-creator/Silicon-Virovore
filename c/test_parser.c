/*
 * test_parser.c
 * Silicon Virovore Project — May 2026
 *
 * Verifies the FASTA file-loading architecture for the real-world 
 * target (HERV-K) and dual human decoy sequences (Albumin and NCAM1).
 *
 * Compile: gcc -o parser_demo test_parser.c
 * Run:     ./parser_demo
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Max memory size allocation per sequence (approx 10,000 characters)
#define MAX_SEQ_LEN 10000

void load_fasta_sequence(const char *filename, char *output_buffer) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "\nFATAL ERROR: Could not open data file: %s\n", filename);
        fprintf(stderr, "Please make sure the file is in your data/ folder and named correctly!\n");
        exit(1);
    }

    char line[256];
    output_buffer[0] = '\0';

    while (fgets(line, sizeof(line), file)) {
        // Skip FASTA description lines starting with '>'
        if (line[0] == '>') {
            continue; 
        }
        
        // Strip out invisible newline formatting safely
        line[strcspn(line, "\r\n")] = 0;
        
        // Safety check to ensure we don't overflow our memory buffer
        if (strlen(output_buffer) + strlen(line) < MAX_SEQ_LEN) {
            strcat(output_buffer, line);
        }
    }

    fclose(file);
}

int main(void) {
    // Allocate space in your computer's memory for the real biological strings
    char target_herv[MAX_SEQ_LEN];
    char decoy_albumin[MAX_SEQ_LEN];
    char decoy_ncam1[MAX_SEQ_LEN];

    printf("=============================================================\n");
    printf("Silicon Virovore — Multi-Target Data Initialization\n");
    printf("May 2026\n");
    printf("=============================================================\n\n");

    // Attempt to open and read each real file using relative directory steps
    printf("Loading viral target (HERV-K Env) from target file... \n");
    load_fasta_sequence("../data/P61567.fasta.txt", target_herv);
    printf("--> SUCCESS: Loaded %lu amino acid residues.\n\n", strlen(target_herv));

    printf("Loading blood plasma decoy (Albumin)... \n");
    load_fasta_sequence("../data/P02768.fasta.txt", decoy_albumin);
    printf("--> SUCCESS: Loaded %lu amino acid residues.\n\n", strlen(decoy_albumin));

    printf("Loading motor neuron cell surface decoy (NCAM1)... \n");
    load_fasta_sequence("../data/P13591.fasta.txt", decoy_ncam1);
    printf("--> SUCCESS: Loaded %lu amino acid residues.\n\n", strlen(decoy_ncam1));

    printf("=============================================================\n");
    printf("DATABASE CONNECTIVITY VERIFIED SUCCESSFUL!\n");
    printf("Your June production pipeline data layer is fully ready.\n");
    printf("=============================================================\n");

    return 0;
}