#include <stdio.h>
#include <string.h>
#include <ctype.h>

// Helper to check standard single-letter amino acid codes
static int is_valid_amino_acid(char c) {
    const char *valid_aa = "ACDEFGHIKLMNPQRSTVWY";
    return (strchr(valid_aa, toupper((unsigned char)c)) != NULL);
}

double c_check_sequence_fitness(const char* sequence) {
    if (sequence == NULL) {
        return -999.0;
    }
    
    int len = strlen(sequence);
    if (len < 9) {
        return -1.0;
    }
    
    // Validate characters to prevent Segfaults
    for (int i = 0; i < len; i++) {
        if (!is_valid_amino_acid(sequence[i])) {
            fprintf(stderr, "[C ERROR] Invalid amino acid '%c' at position %d\n", sequence[i], i);
            return -1.0; // Return invalid score safely without crashing C
        }
    }
    
    // Wire your internal C scoring logic safely here
    return 1.0; 
}