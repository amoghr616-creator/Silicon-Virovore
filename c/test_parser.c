// test_parser.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "virovore.h"

Variant init_variant(size_t initial_capacity) {
    Variant v;
    v.sequence = (char *)malloc(initial_capacity * sizeof(char));
    if (!v.sequence) {
        fprintf(stderr, "Heap allocation failed.\n");
        exit(1);
    }
    v.length = 0;
    v.capacity = initial_capacity;
    v.sequence[0] = '\0';
    return v;
}

void append_sequence(Variant *v, const char *chunk) {
    size_t chunk_len = strlen(chunk);
    
    // Capacity doubling loop logic to guarantee zero fragmentation overhead
    if (v->length + chunk_len + 1 > v->capacity) {
        v->capacity = (v->capacity * 2) + chunk_len; 
        v->sequence = (char *)realloc(v->sequence, v->capacity * sizeof(char));
        if (!v->sequence) {
            fprintf(stderr, "Realloc fragmentation failure.\n");
            exit(1);
        }
    }
    
    strcat(v->sequence, chunk);
    v->length += chunk_len;
}