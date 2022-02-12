#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* inspectfa.c - Version 1 */


int default_samp_rate = 16;

int testing_lite = 0;


int main(int argc, char *argv[]) {

    /* INPUT VARS */
    char *index_path;
    int sample_rate;
    char *output_path;

    /* OUTPUT VARS */
    char *bwt;

    int num_a;
    int num_c;
    int num_g;
    int num_t;

    int *tally_a;
    int *tally_c;
    int *tally_g;
    int *tally_t;

    /* OTHER VARS */
    FILE *index_file;
    FILE *output_file;

    int genome_len = -1;

    char *text_buffer;
    char *int_buffer;

    char *genome;

    int *suffix_arr;


    int x;


    /* GET INPUTS */
    if (argc == 4) {
        index_path = argv[1];
        sample_rate = atoi(argv[2]);
        output_path = argv[3];
    } else {
        index_path = "default_FM";
        sample_rate = default_samp_rate;
        output_path = "default_inspect";
    }

    /* READ FROM BINARY FILE - GET SUFFIX ARRAY */
    index_file = fopen(index_path, "rb");

    genome_len = 0;
    /* Read the genome_len */
    fread(&genome_len, sizeof(genome_len), 1, index_file);

    if(testing_lite) {
        printf("genome length = %d\n", genome_len);
    }

    /* Allocate space for genome string and suffix_arr */
    genome = (char*) malloc(genome_len + 4);
    suffix_arr = (int*) malloc(genome_len * sizeof(int));

    /* read the genome */
    fread(genome, sizeof(char), genome_len, index_file);

    /* Add null terminator to genome - allows for string processing */
    genome[genome_len] = '\0';

    /* Confirm string length */
    if(testing_lite) {
        printf("genome string length = %d\n", (int)strlen(genome));
    }
    
    //printf("%s\n", genome);

    /* get suffix_arr */
    fread(suffix_arr, sizeof(int), genome_len, index_file);

    /* Print suffix_arr to check */
    /*
    printf("Suffix array (by 100s):\n");
    for (int i = 0 ; i < genome_len ; i += 100){
        printf("%d\n", suffix_arr[i]);
    }
    */

    /* READ ADDITIONAL INFORMATION FOR FM INDEX */
    bwt = malloc((genome_len + 4) * sizeof(char));
    tally_a = malloc( (genome_len + 1) * sizeof(int));
    tally_c = malloc( (genome_len + 1) * sizeof(int));
    tally_g = malloc( (genome_len + 1) * sizeof(int));
    tally_t = malloc( (genome_len + 1) * sizeof(int));

    fread(bwt, sizeof(char), genome_len, index_file);

    fread(&num_a, sizeof(int), 1, index_file);
    fread(&num_c, sizeof(int), 1, index_file);
    fread(&num_g, sizeof(int), 1, index_file);
    fread(&num_t, sizeof(int), 1, index_file);

    fread(tally_a, sizeof(int), genome_len, index_file);
    fread(tally_c, sizeof(int), genome_len, index_file);
    fread(tally_g, sizeof(int), genome_len, index_file);
    fread(tally_t, sizeof(int), genome_len, index_file);



    /* WRITE SUFFIX ARRAY VALUES TO OUTPUT FILE */
    output_file = fopen(output_path, "w");

    text_buffer = malloc ( 100 * sizeof(char) );

    sprintf(text_buffer, "1\t%d\t%d\t%d\t%d\t\n", num_a, num_c, num_g, num_t);

    if(testing_lite) {
        printf("%s", text_buffer);
    }
    fputs(text_buffer, output_file);

    free(text_buffer);

    fputs(bwt, output_file);
    fputs("\n", output_file);

    /* Write selected entries from tally table */
    int_buffer = malloc(16 * sizeof(char));

    /* TALLY SPOT CHECK A */
    x = 0;
    while ( x <= genome_len / sample_rate) {

        sprintf(int_buffer, "%d\t", tally_a[x * sample_rate]);

        fputs(int_buffer, output_file);

        //printf("%s\n", int_buffer);

        x++;
    }
    fputs("\n", output_file);

    /* TALLY SPOT CHECK C */
    x = 0;
    while ( x <= genome_len / sample_rate) {

        sprintf(int_buffer, "%d\t", tally_c[x * sample_rate]);

        fputs(int_buffer, output_file);

        //printf("%s\n", int_buffer);

        x++;
    }
    fputs("\n", output_file);

    /* TALLY SPOT CHECK G */
    x = 0;
    while ( x <= genome_len / sample_rate) {

        sprintf(int_buffer, "%d\t", tally_g[x * sample_rate]);

        fputs(int_buffer, output_file);

        //printf("%s\n", int_buffer);

        x++;
    }
    fputs("\n", output_file);

    /* TALLY SPOT CHECK T */
    x = 0;
    while ( x <= genome_len / sample_rate) {

        sprintf(int_buffer, "%d\t", tally_t[x * sample_rate]);

        fputs(int_buffer, output_file);

        //printf("%s\n", int_buffer);

        x++;
    }
    fputs("\n", output_file);

    free(int_buffer);

    //printf("%s\n", text_buffer);
    //fputs(text_buffer, output_file);

    fclose(output_file);

    /* Free pointers used */
    free(suffix_arr);
    free(genome);

    free(bwt);
    free(tally_a);
    free(tally_c);
    free(tally_g);
    free(tally_t);

    printf("Done running inspectfm\n^-^\n");

    return 0;
}