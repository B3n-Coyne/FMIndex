#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* buildfm.c - Version 1 */

int testing_lite = 0;

int testing_fm = 0;

int testing_fm_full = 0;

/* Changes which imput file will be used as the default when no arguments are given */
int default_input = 1;
// 1 = ecoli_tiny.txt (small genome)
// 0 = ecoli.txt (full genome)

int leq2(int a1, int a2, int b1, int b2) {
    return(a1 < b1 || a1 == b1 && a2 <= b2);
}

int leq3(int a1, int a2, int a3, int b1, int b2, int b3) {
    return(a1 < b1 || a1 == b1 && leq2(a2,a3, b2,b3));
}

void radix_pass(int* a, int* b, int* r, int n, int K) {
    int* c = malloc((K + 1) * sizeof(int)); // counter array
    for (int i = 0; i <= K; i++) {
        c[i] = 0; // reset counters
    }

    for (int i = 0; i < n; i++) {
        c[r[a[i]]]++; // count occurrences
    }

    for (int i = 0, sum = 0; i <= K; i++) { // exclusive prefix sums
        int t = c[i];
        c[i] = sum;
        sum += t;
    }
    for (int i = 0; i < n; i++) { // sort
        b[c[r[a[i]]]++] = a[i];
    }

    free(c);
}

/* Recursive DC3 algorithm */
int* DC3_construct( int *s, int *SA, int n, int K) {
    int n0=(n+2)/3, n1=(n+1)/3, n2=n/3;
    int n02=n0+n2;

    int* s12 = malloc((n02+3) * sizeof(int) ); //new int[n02 + 3];
    s12[n02] = 0;
    s12[n02+1] = 0;
    s12[n02+2] = 0;

    int* SA12 = malloc((n02+3) * sizeof(int) ); //new int[n02 + 3];
    SA12[n02] = 0;
    SA12[n02+1] = 0;
    SA12[n02+2] = 0;

    int* s0 = malloc(n0 * sizeof(int)); //new int[n0];
    int* SA0 = malloc(n0 * sizeof(int)); //new int[n0];

    // generate positions of mod 1 and mod 2 suffixes
    // the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
    for (int i=0, j=0; i < n+(n0-n1); i++) {
        if (i%3 != 0) {
            s12[j++] = i;
        }
    }

    // lsb radix sort the mod 1 and mod 2 triples
    radix_pass(s12 , SA12, s+2, n02, K);
    radix_pass(SA12, s12 , s+1, n02, K);
    radix_pass(s12 , SA12, s , n02, K);

    // find lexicographic names of triples
    int name = 0, c0 = -1, c1 = -1, c2 = -1;
    for (int i = 0; i < n02; i++) {
        if (s[SA12[i]] != c0 || s[SA12[i]+1] != c1 || s[SA12[i]+2] != c2) {
            name++;
            c0 = s[SA12[i]];
            c1 = s[SA12[i]+1];
            c2 = s[SA12[i]+2];
        }

        if (SA12[i] % 3 == 1) { // left half
            s12[SA12[i]/3] = name;
        } else { // right half
            s12[SA12[i]/3 + n0] = name;
        } 
    }

    // recurse if names are not yet unique
    if (name < n02) {
        DC3_construct(s12, SA12, n02, name);

        // store unique names in s12 using the suffix array
        for (int i = 0; i < n02; i++) {
            s12[SA12[i]] = i + 1;
        }

    } else { // generate the suffix array of s12 directly
        for (int i = 0; i < n02; i++) {
            SA12[s12[i] - 1] = i;
        }
    }

    // stably sort the mod 0 suffixes from SA12 by their first character
    for (int i=0, j=0; i < n02; i++) {
        if (SA12[i] < n0) {
            s0[j++] = 3*SA12[i];
        }
    }
    radix_pass(s0, SA0, s, n0, K);

    // merge sorted SA0 suffixes and sorted SA12 suffixes
    for (int p=0, t=n0-n1, k=0; k < n; k++) {
        #define GetI() (SA12[t] < n0 ? SA12[t]*3+1: (SA12[t] - n0) * 3 + 2)
        int i = GetI(); // pos of current offset 12 suffix
        int j = SA0[p]; // pos of current offset 0 suffix

        // different compares for mod 1 and mod 2 suffixes
        if (SA12[t] < n0 ?  leq2(s[i], s12[SA12[t] + n0], s[j], s12[j/3]) :
                            leq3(s[i],s[i+1],s12[SA12[t]-n0+1], s[j],s[j+1],s12[j/3+n0]) ) {// suffix from SA12 is smaller
            SA[k] = i; t++;
            if (t == n02) { // done --- only SA0 suffixes left
                for (k++; p < n0; p++, k++) SA[k] = SA0[p];
            }
        } else {// suffix from SA0 is smaller
            SA[k] = j; p++;
            if (p == n0) { // done --- only SA12 suffixes left
                for (k++; t < n02; t++, k++) SA[k] = GetI();
            }
        }
    }
    //delete [] s12; delete [] SA12; delete [] SA0; delete [] s0;
    free(s12);
    free(SA12);
    free(SA0);
    free(s0);

    return SA;
}

/* Outer function to construct the suffix array - Calls recursive DC3 function */
int* SAconstruct(char *ref_text, int genome_len) {

    int *suffix_arr;

    int *ranks;

    int *samp_suffixes;
    int *non_samp_suffixes;

    char *padded_genome;
    int *int_genome;

    int padded_len;

    int samp_i;
    int non_samp_i;

    padded_genome = malloc(  (genome_len + 3) * sizeof(int) );

    padded_len = genome_len;

    strcpy(padded_genome, ref_text);

    if (strlen(padded_genome) % 3 != 0 ) {
        strcat(padded_genome, "$");
        padded_len++;
    }
    /* length of padded_genome should now be a multiple of 3 */
    //padded_len = strlen(padded_genome);

    /* Map padded_genome to ints in int_genome for easier sorting */
    int_genome = malloc( padded_len * sizeof(int) );
    for ( int i = 0 ; i < padded_len ; i++ ) {
        if(padded_genome[i] == 'A') {
            int_genome[i] = 1;
        } else if (padded_genome[i] == 'C') {
            int_genome[i] = 2;
        } else if (padded_genome[i] == 'G') {
            int_genome[i] = 3;
        } else if (padded_genome[i] == 'T') {
            int_genome[i] = 4;
        } else { // Teminal char: $
            int_genome[i] = 0;
        }
    }

    /* Want suffix array the length of the padded array */
    suffix_arr = malloc( (padded_len + 1) * sizeof(int) );
    if(suffix_arr == NULL) {
        printf("Could not allocate space for suffix_arr\n");
        return NULL;
    }

    /* Initialize suffix array */
    for ( int i = 0 ; i < padded_len ; i++ ) {
        suffix_arr[i] = -1;
    }
    
    if(testing_lite) {
        printf("Beginning suffix array construction\n");
    }   

    /* Not really sure what K (last parameter) should be here... - Setting as genome_len, try padded_len? */
    suffix_arr = DC3_construct( int_genome, suffix_arr, padded_len, padded_len);

    if(testing_lite) {
        printf("Done running DC3 algorithm\n");
    }

    /* Want to trim 1 char from front of suffix array for every character added in padding */
    while (padded_len > genome_len) {
        for ( int i = 1 ; i < padded_len ; i++ ) {
            suffix_arr[i-1] = suffix_arr[i];
        }
        padded_len--;
    }

    free(padded_genome);
    free(int_genome);

    return suffix_arr;

}

/* Main function to process string/genome */

int main(int argc, char *argv[]) {

    /* INPUT VARS */
    char *reference_path;
    char *output_path;

    /* OTHER VARS */
    char *genome;

    int *suffix_arr;

    char *text_buffer;

    int genome_len = -1;

    int i;

    FILE *reference_file;
    FILE *output_file;

    /* GET INPUTS */
    if (argc == 3) {
        reference_path = argv[1];
        output_path = argv[2];
    } else {
        if (default_input) {
            reference_path = "ecoli_tiny.txt";
        } else {
            reference_path = "ecoli.txt";
        }
        output_path = "default_FM";
    }

    /* GET THE GENOME/TEXT */

    reference_file = fopen(reference_path, "r");

    /* Max genome size of 100,000,000 */
    genome = malloc(100000000 * sizeof(char));
    if (genome == NULL) {
        printf("Could not allocate space for genome");
        return 0;
    }
    genome[0] = '\0';

    text_buffer = malloc(10000 * sizeof(char));
    if (text_buffer == NULL) {
        printf("Could not allocate space for text_buffer (1)");
        return 0;
    }

    /* START READING FILE */
    text_buffer = fgets(text_buffer, 10000, reference_file);

    if(testing_lite) {
        printf("Header: %s\n", text_buffer);
    }

    while (!feof(reference_file)) {
        /* Read next line in */
        fscanf(reference_file, "%s\n", text_buffer);

        strcat(genome, text_buffer);
    }

    free(text_buffer);

    /*
    printf("%s\n", genome);
    */

    /* Using '$' as terinal character */
    strcat(genome, "$");

    genome_len = strlen(genome);

    if(testing_lite) {
        printf("genome_len = %d\n", genome_len);
    }


    /* SORT/CREATE SUFFIX ARRAY */
    suffix_arr = SAconstruct(genome, genome_len);

    if (suffix_arr == NULL) {
        printf("Something went wrong constructing the suffix array");
        return 0;
    }

    /*
    printf("Suffix array (by 100s):\n");
    for (int i = 0 ; i < genome_len ; i += 100){
        printf("%d\n", suffix_arr[i]);
    }
    */

    /****************/
    /* GET FM INDEX */
    /****************/
    /* WANT:
     *      BWT of the genome
     *      4 ints to count occurences of each char at beginning of FM index line (first column)
     *      Tally table for occ queries
     */

    if (testing_fm) {
        printf("Beginning FM Index construction\n");
    }

    char *bwt;
    int num_a = 0, num_c = 0, num_g = 0, num_t = 0;

    int *tally_a;
    int *tally_c;
    int *tally_g;
    int *tally_t;

    int temp;

    /* Initialize values for fm index */

    bwt = malloc(genome_len * sizeof(char));
    if (bwt == NULL) {
        printf("Could not allocate space for the BWT\n");
        return 0;
    }

    tally_a = malloc(genome_len * sizeof(int));
    tally_c = malloc(genome_len * sizeof(int));
    tally_g = malloc(genome_len * sizeof(int));
    tally_t = malloc(genome_len * sizeof(int));

    for (int i = 0 ; i < genome_len ; i++ ) {
        tally_a[i] = 0;
        tally_c[i] = 0;
        tally_g[i] = 0;
        tally_t[i] = 0;
    }

    /* Create the BWT from the suffix array (and the first column of the matrix) */

    bwt[0] = genome[genome_len - 2];

    /* Get first value for tally array */
    if (bwt[0] == 'a' || bwt[0] == 'A') {
        tally_a[0] = 1;
    } else if (bwt[0] == 'c' || bwt[0] == 'C') {
        tally_c[0] = 1;
    } else if (bwt[0] == 'g' || bwt[0] == 'G') {
        tally_g[0] = 1;
    } else if (bwt[0] == 't' || bwt[0] == 'T') {
        tally_t[0] = 1;
    }

    for (int i = 1 ; i < genome_len ; i++) {

        /* Create the BWT */
        temp = suffix_arr[i] - 1;

        if (temp < 0) {
            temp += genome_len;
        }

        bwt[i] = genome[temp];

        /* Count for the first column */
        if ( genome[suffix_arr[i]] == 'a' || genome[suffix_arr[i]] == 'A') {
            num_a++;
        } else if ( genome[suffix_arr[i]] == 'c' || genome[suffix_arr[i]] == 'C' ) {
            num_c++;
        } else if ( genome[suffix_arr[i]] == 'g' || genome[suffix_arr[i]] == 'G' ) {
            num_g++;
        } else if ( genome[suffix_arr[i]] == 't' || genome[suffix_arr[i]] == 'T' ) {
            num_t++;
        }

        /* Update tally array */
        tally_a[i] = tally_a[i - 1];
        tally_c[i] = tally_c[i - 1];
        tally_g[i] = tally_g[i - 1];
        tally_t[i] = tally_t[i - 1];

        if (bwt[i] == 'a' || bwt[i] == 'A') {
            tally_a[i] = tally_a[i] + 1;
        } else if (bwt[i] == 'c' || bwt[i] == 'C') {
            tally_c[i] = tally_c[i] + 1;
        } else if (bwt[i] == 'g' || bwt[i] == 'G') {
            tally_g[i] = tally_g[i] + 1;
        } else if (bwt[i] == 't' || bwt[i] == 'T') {
            tally_t[i] = tally_t[i] + 1;
        }

    }

    if(testing_fm_full) {
        for (int i = 0 ; i < genome_len ; i++ ) {
            printf("%d:\t%c\t%d\t%d\t%d\t%d\t\n", i, bwt[i], tally_a[i], tally_c[i], tally_g[i], tally_t[i]);
        }

        //printf("%s\n", genome);

        printf("%d:\t%c\t%d\t%d\t%d\t%d\t\n", 0, bwt[0], tally_a[0], tally_c[0], tally_g[0], tally_t[0]);
    }

    if (testing_fm) {
        printf("Done constructing FM Index\n");
    }


    /* WRITE INFO TO OUTPUT FILE */
    output_file = fopen(output_path, "wb");

    /* Need in output file: genome, genome_len, sorted suffix_arr */
    
    /* Write length of genome */
    fwrite(&genome_len, sizeof(genome_len), 1, output_file);

    /* Write genome (no newline '\0') */
    fwrite(genome, sizeof(char), genome_len, output_file);

    /* Write suffix_arr */
    fwrite(suffix_arr, sizeof(int), genome_len, output_file);

    /******************************/
    /* WRITE FM INDEXING INFO HERE*/
    /******************************/
    /* WANT:
     *      BWT of the genome
     *      4 ints to count occurences of each char at beginning of FM index line (first column)
     *      Tally table for occ queries
     */

    fwrite(bwt, sizeof(char), genome_len, output_file);
    fwrite(&num_a, sizeof(int), 1, output_file);
    fwrite(&num_c, sizeof(int), 1, output_file);
    fwrite(&num_g, sizeof(int), 1, output_file);
    fwrite(&num_t, sizeof(int), 1, output_file);

    fwrite(tally_a, sizeof(int), genome_len, output_file);
    fwrite(tally_c, sizeof(int), genome_len, output_file);
    fwrite(tally_g, sizeof(int), genome_len, output_file);
    fwrite(tally_t, sizeof(int), genome_len, output_file);

    fclose(output_file);

    free(suffix_arr);

    free(genome);

    free(bwt);
    free(tally_a);
    free(tally_g);
    free(tally_c);
    free(tally_t);

    printf("Done running buildfm\n^-^\n");

    return 0;
}