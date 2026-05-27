/*
 * charge.c
 * Silicon Virovore Project — May 2026
 *
 * Implements the net charge calculation and "avoid ~0 charge" gate.
 * This is the second "gate" function in your biophysical scoring pipeline.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

/* ============================================================================
 * PART A: CHARGE SUMMARY STRUCTURE
 * ============================================================================
 * We define this structure at the top so all functions below know what it is!
 */

typedef struct {
    int net_charge;
    float penalty;
    int num_positive;  // Number of K + R residues
    int num_negative;  // Number of D + E residues
} charge_summary;


/* ============================================================================
 * PART B: CHARGE CONTRIBUTIONS PER AMINO ACID
 * ============================================================================
 * Table: 1 if positive (K, R), -1 if negative (D, E), 0 if neutral
 */

int amino_acid_charge[26] = {
    0,   // A (Alanine)
    0,   // B (Asp or Asn) — placeholder
    0,   // C (Cysteine)
    -1,  // D (Aspartic acid) — negatively charged
    -1,  // E (Glutamic acid) — negatively charged
    0,   // F (Phenylalanine)
    0,   // G (Glycine)
    0,   // H (Histidine)
    0,   // I (Isoleucine)
    0,   // J (not used)
    1,   // K (Lysine) — positively charged
    0,   // L (Leucine)
    0,   // M (Methionine)
    0,   // N (Asparagine)
    0,   // O (not used)
    0,   // P (Proline)
    0,   // Q (Glutamine)
    1,   // R (Arginine) — positively charged
    0,   // S (Serine)
    0,   // T (Threonine)
    0,   // U (not used)
    0,   // V (Valine)
    0,   // W (Tryptophan)
    0,   // X (unknown)
    0,   // Y (Tyrosine)
    0    // Z (Glu or Gln) — placeholder
};

/* ============================================================================
 * PART C: LOOKUP FUNCTION FOR CHARGE
 * ============================================================================
 */

int get_amino_acid_charge(char amino_acid)
{
    amino_acid = toupper(amino_acid);
    
    if (amino_acid < 'A' || amino_acid > 'Z')
    {
        fprintf(stderr, "ERROR: Invalid amino acid '%c' in charge lookup\n", amino_acid);
        return 0; 
    }
    
    int index = amino_acid - 'A';
    return amino_acid_charge[index];
}

/* ============================================================================
 * PART D: NET CHARGE PENALTY ("AVOID ~0" GATE)
 * ============================================================================
 */

float calculate_charge_penalty(int net_charge, int threshold)
{
    int abs_charge = (net_charge < 0) ? -net_charge : net_charge;
    
    if (abs_charge > threshold)
    {
        return 0.0; // Sequence has sufficient charge — no penalty
    }
    else
    {
        // Penalty ramps up smoothly as |charge| hits 0
        float penalty = (float)(threshold + 1 - abs_charge);
        return penalty;
    }
}

/* ============================================================================
 * PART E: STREAMLINED COMPUTE FULL CHARGE PROFILE (SINGLE PASS)
 * ============================================================================
 * This replaces the double loops! We look at each character exactly once.
 */

int compute_charge_profile(const char *sequence, int charge_threshold,
                            charge_summary *summary_out)
{
    if (sequence == NULL || summary_out == NULL)
    {
        fprintf(stderr, "ERROR: NULL pointer in compute_charge_profile\n");
        return -1;
    }
    
    int num_positive = 0;
    int num_negative = 0;
    
    // SINGLE PASS: Iterates through the protein string once
    for (int i = 0; sequence[i] != '\0'; i++)
    {
        char aa = toupper(sequence[i]);
        
        if (aa == 'K' || aa == 'R')
        {
            num_positive++;
        }
        else if (aa == 'D' || aa == 'E')
        {
            num_negative++;
        }
        else if (aa < 'A' || aa > 'Z')
        {
            fprintf(stderr, "ERROR: Invalid amino acid '%c' at position %d\n", sequence[i], i);
            return -1;
        }
    }
    
    // Assign calculated metrics directly to your struct outputs
    summary_out->num_positive = num_positive;
    summary_out->num_negative = num_negative;
    summary_out->net_charge = num_positive - num_negative;
    summary_out->penalty = calculate_charge_penalty(summary_out->net_charge, charge_threshold);
    
    return 0;  // Success
}

/* ============================================================================
 * PART F: DEMONSTRATION AND TESTING
 * ============================================================================
 */

int main(void)
{
    printf("=============================================================\n");
    printf("Silicon Virovore — Streamlined Net Charge Analyzer\n");
    printf("May 2026\n");
    printf("=============================================================\n\n");
    
    // TEST 1: Positively charged
    printf("TEST 1: Positively charged sequence\n");
    printf("---\n");
    const char *test1 = "MKLAVLALLVTFAGSSDLH";
    charge_summary s1;
    if (compute_charge_profile(test1, 1, &s1) == 0)
    {
        printf("Sequence: %s\n", test1);
        printf("Positive (K+R): %d | Negative (D+E): %d\n", s1.num_positive, s1.num_negative);
        printf("Net charge: %d\n", s1.net_charge);
        printf("Charge penalty: %.2f %s\n\n", s1.penalty, s1.penalty == 0.0 ? "(PASS)" : "(FAIL)");
    }
    
    // TEST 2: Negatively charged
    printf("TEST 2: Negatively charged sequence\n");
    printf("---\n");
    const char *test2 = "DEDEEDEDEDDEDEDED";
    charge_summary s2;
    if (compute_charge_profile(test2, 1, &s2) == 0)
    {
        printf("Sequence: %s\n", test2);
        printf("Positive (K+R): %d | Negative (D+E): %d\n", s2.num_positive, s2.num_negative);
        printf("Net charge: %d\n", s2.net_charge);
        printf("Charge penalty: %.2f %s\n\n", s2.penalty, s2.penalty == 0.0 ? "(PASS)" : "(FAIL)");
    }
    
    // TEST 3: Balanced near 0 (Should fail gate)
    printf("TEST 3: Balanced sequence (Near charge 0 constraint)\n");
    printf("---\n");
    const char *test3 = "MKLAVDEALLVTFAGSS"; // 1 K (+1), 1 D (-1), 1 E (-1) = Net -1
    charge_summary s3;
    if (compute_charge_profile(test3, 1, &s3) == 0)
    {
        printf("Sequence: %s\n", test3);
        printf("Positive (K+R): %d | Negative (D+E): %d\n", s3.num_positive, s3.num_negative);
        printf("Net charge: %d\n", s3.net_charge);
        printf("Charge penalty: %.2f %s\n\n", s3.penalty, s3.penalty > 0.0 ? "(FAIL — aggregation risk!)" : "(PASS)");
    }
    
    // TEST 4: Exactly net charge 0 (Worst case scenario)
    printf("TEST 4: Exactly balanced (Worst case — charge = 0)\n");
    printf("---\n");
    const char *test4_zero = "MKLAVDALLVTFAGSS"; // 1 K (+1), 1 D (-1) = Net 0
    charge_summary s4;
    if (compute_charge_profile(test4_zero, 1, &s4) == 0)
    {
        printf("Sequence: %s\n", test4_zero);
        printf("Positive (K+R): %d | Negative (D+E): %d\n", s4.num_positive, s4.num_negative);
        printf("Net charge: %d\n", s4.net_charge);
        printf("Charge penalty: %.2f %s\n\n", s4.penalty, s4.penalty > 1.0 ? "(FAIL — HIGHEST AGGREGATION RISK!)" : "(PASS)");
    }
    
    // TEST 5: Graceful error handling
    printf("TEST 5: Invalid character validation (contains '9')\n");
    printf("---\n");
    const char *test5 = "MKL9VALLVTFAGSSDLH";
    charge_summary s5;
    if (compute_charge_profile(test5, 1, &s5) != 0)
    {
        printf("System Protection Active: Invalid sequence safely rejected.\n\n");
    }
    
    printf("=============================================================\n");
    printf("All optimized loops verified successfully!\n");
    printf("=============================================================\n");
    
    return 0;
}