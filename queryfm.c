#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* queryfm.c - Version 1 */

/* MODIFY FOR FM INDEXING - WILL NEED TO BE SIGNIFICANTLY MODIFIED*/

/* DIAGNOSTIC VARS */
// Used to confirm proper program flow - prints after major computations
int testing_light = 0;
// Used to confirm data - prints all major data calculated
// WARNING: Will print full data - avoid with large datasets!
int testing_full  = 0;

// Used to select default search method
int default_search = 0;
// 0 = partial
// 1 = complete


/* Main function to process queries of a suffix array */
int main(int argc, char *argv[]) {

    

    /* INPUT VARS */
    char *index_path;
    char *queries_path;
    //char *query_mode;
    char *output_path;

    /* OTHER VARS */

    char *genome;

    char *text_buffer;

    int *hit_locations;

    int *suffix_arr;

    int genome_len;

    int num_hits;

    FILE *index_file;
    FILE *output_file;
    FILE *query_file;

    int hit_header;

    char **header_arr;
    char **query_arr;

    int num_queries = 0;

    char *num_buffer;

    int comp_search;


    /* FM Indexing vars */
    char *bwt;
    int num_a, num_c, num_g, num_t;
    int *tally_a;
    int *tally_c;
    int *tally_g;
    int *tally_t;


    if (testing_light) {
        printf ("\nBEGINNING RUN OF querysa\n\n");
    }

    /* GET INPUTS */
    if (argc == 5) {
        index_path = argv[1];
        queries_path = argv[2];
        //query_mode = argv[3];
        output_path = argv[4];


        /* "complete" or "partial" */
        if (strcmp(argv[3], "partial") == 0) {
            comp_search = 0;
        } else {
            comp_search = 1;
        }
    } else {
        index_path = "default_FM";
        queries_path = "reads.txt";
        //query_mode = "naive";
        output_path = "default_queries";

        comp_search = default_search;
    }

    /* READ FROM BINARY FILE - GET SUFFIX ARRAY */
    index_file = fopen(index_path, "rb");

    genome_len = 0;
    /* Read the genome_len */
    fread(&genome_len, sizeof(genome_len), 1, index_file);
    if (testing_light) {
        printf("genome_len = %d\n", genome_len);
    }

    /* Allocate space for genome string and suffix_arr */
    genome = (char*) malloc(genome_len + 1);
    suffix_arr = (int*) malloc(genome_len * sizeof(int));

    /* read the genome */
    fread(genome, sizeof(char), genome_len, index_file);

    /* Add null terminator to genome - allows for string processing */
    genome[genome_len] = '\0';

    /* Confirm string length */
    if (testing_light) {
        printf("strlen(genome) = %d\n", (int)strlen(genome));
    }
    
    if (testing_full) {
        printf("Genome:\n%s\n", genome);
    }

    /* get suffix_arr */
    fread(suffix_arr, sizeof(int), genome_len, index_file);

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


    output_file = fopen(output_path, "w");
    query_file = fopen(queries_path, "r");


    text_buffer = malloc(100000 * sizeof(char) );
    header_arr = malloc(10001 * sizeof(char*) );
    query_arr = malloc(10000 * sizeof(char*));

    /* Read first header */
    header_arr[0] = malloc(100 * sizeof(char));
    fscanf(query_file, ">%s\n", header_arr[0]);

    num_queries = 0;

    /* GET ARRAY OF HEADERS AND QUERIES */
    while ( !feof(query_file) && text_buffer != NULL) {

        if (testing_full) {
            printf("Reading Query %d\n", num_queries);
        }

        query_arr[num_queries] = malloc(100000 * sizeof(char));
        query_arr[num_queries][0] = '\0';

        header_arr[num_queries + 1] = malloc(100 * sizeof(char));

        hit_header = 0;


        while (!hit_header) {
            fscanf(query_file, "%s\n", text_buffer);

            if (sscanf(text_buffer, ">%s", header_arr[num_queries + 1]) == 1 ) {
                hit_header = 1;
            } else {
                strcat(query_arr[num_queries], text_buffer);
            }

            if (feof(query_file)) {
                hit_header = 1;
            }
        }

        num_queries++;
    }
        
    if (testing_light) {
        printf("num_queries = %d\n", num_queries);
    }



    /* FOR EACH QUERY */
    for (int i = 0 ; i < num_queries ; i++) {

        if (testing_full) {
            printf("%s\n%s\n", header_arr[i], query_arr[i]);
        }

        /* WRITE HEADER TO text_buffer */
        strcpy(text_buffer, header_arr[i]);
        strcat(text_buffer, "\t");

        fputs(text_buffer, output_file);
        
        num_hits = 0;

        hit_locations = malloc(genome_len * sizeof(int));
        /* Assuming for hit_locations that there aren't going to be more than 100,000 occurences of a substring in the genome */

        //int init_hit = -1;

        /*
        int right = genome_len - 1;
        int left = 0;
        int center = genome_len / 2;
        */

        int query_len = strlen(query_arr[i]);

        

        /* SEARCH HERE - FM INDEX QUERY */
        
        /* if comp_search - do full search thru index, saving each full match that is found */

        /* if !comp_search - do full search thru index, if find a new partial match that is longer than previous, save it as the new longest, reset the found hits
         * If find a match that is shorter than the previous longest, disregard,
         * If find a match that is the same length as the previous longest, save it as an additional hit
         */

        /* GENERAL ALGORITHM: Going thru string in reverse order
         * Use first column to find range which last char in query occurs
         * Look at what character is next in query, go to that occ tally table
         * In tally table, find the tally at beginning and end - use difference, offset by first tally, to find range in first column
         * Repeat for furthur characters, using updated range
         * If the query has been completed, all indicies in the range are complete hits
         * If the query cannot be completed, all indicies in the previous range are partial hits
         */

        num_hits = 0;

        int right = 0;
        int left = 0;

        int hit_len;

        //char tally_list = 'a';

        int curr_char = query_len - 1;

        /* Get starting range */
        if (query_arr[i][curr_char] == 'a' || query_arr[i][curr_char] == 'A') {
            right = 1; //Skipping the '$' in the first row
            left = right + num_a;
        } else if (query_arr[i][curr_char] == 'c' || query_arr[i][curr_char] == 'C') {
            right = 1 + num_a;
            left = right + num_c;
        } else if (query_arr[i][curr_char] == 'g' || query_arr[i][curr_char] == 'G') {
            right = 1 + num_a + num_c;
            left = right + num_g;
        } else if (query_arr[i][curr_char] == 't' || query_arr[i][curr_char] == 'T') {
            right = 1 + num_a + num_c + num_g;
            left = right + num_t;
        }

        /* num_hits will be the size of the final range */

        int found_hits = 0;

        hit_len = 1;

        curr_char--;

        while (!found_hits) {

            /* Use current range to check tally table for next character */
            /* Use values found in tally table to get new range */
            int smol_tally = -1;
            int big_tally = -1;

            if (query_arr[i][curr_char] == 'a' || query_arr[i][curr_char] == 'A') {
                smol_tally = tally_a[right - 1];
                big_tally = tally_a[left - 1];
            } else if (query_arr[i][curr_char] == 'c' || query_arr[i][curr_char] == 'C') {
                smol_tally = tally_c[right - 1];
                big_tally = tally_c[left - 1];
            } else if (query_arr[i][curr_char] == 'g' || query_arr[i][curr_char] == 'G') {
                smol_tally = tally_g[right - 1];
                big_tally = tally_g[left - 1];
            } else if (query_arr[i][curr_char] == 't' || query_arr[i][curr_char] == 'T') {
                smol_tally = tally_t[right - 1];
                big_tally = tally_t[left - 1];
            }

            /* If the tally doesn't change in the range */
            if (smol_tally == big_tally) {
                if (!comp_search) {
                    /* Doing a partial search, want to keep the old range */
                    found_hits = 1;
                } else { // Doing a complete search, and didn't find the full query string
                    found_hits = 1;
                    right = -1;
                    left = -1;
                    hit_len = 0;
                }

                break;
            }

            
            hit_len++;

            if (query_arr[i][curr_char] == 'a' || query_arr[i][curr_char] == 'A') {
                right = smol_tally + 1;
                left = big_tally + 1;
            } else if (query_arr[i][curr_char] == 'c' || query_arr[i][curr_char] == 'C') {
                right = smol_tally + 1 + num_a;
                left = big_tally + 1 + num_a;
            } else if (query_arr[i][curr_char] == 'g' || query_arr[i][curr_char] == 'G') {
                right = smol_tally + 1 + num_a + num_c;
                left = big_tally + 1 + num_a + num_c;
            } else if (query_arr[i][curr_char] == 't' || query_arr[i][curr_char] == 'T') {
                right = smol_tally + 1 + num_a + num_c + num_g;
                left = big_tally + 1 + num_a + num_c + num_g;
            }

            curr_char--;

            /* If reached the first character in the query (starting from last) */
            if (curr_char < 0) {
                /* Whatever the current range is, is the final range for complete match */
                found_hits = 1;
                //break;
            }
        }

        /* Use the suffix array values in the found range as the hits */
        num_hits = 0;
        for (int x = right ; x < left ; x++) {
            hit_locations[num_hits] = suffix_arr[x];
            num_hits++;
        }


        /* ADD INFO TO text_buffer */
        num_buffer = malloc(100 * sizeof(char));

        /* Add hit_len */
        sprintf(num_buffer, "%d\t", hit_len);
        //strcat(text_buffer, num_buffer);
        fputs(num_buffer, output_file);

        /* Add num_hits */
        sprintf(num_buffer, "%d\t", num_hits);
        //strcat(text_buffer, num_buffer);
        fputs(num_buffer, output_file);

        /* Add each hit */
        for (int j = 0 ; j < num_hits ; j++) {
            sprintf(num_buffer, "%d\t", hit_locations[j]);

            fputs(num_buffer, output_file);
        }

        fputs("\n", output_file);

        //strcat(text_buffer, "\n");

        free(num_buffer);
        free(hit_locations);

        /* WRITE text_buffer TO output_file */
        //fputs(text_buffer, output_file);
    }

    /* FREE DYNAMICALLY ALLOCATED VARS */

    for (int i = 0 ; i < num_queries ; i++ ){
        free(query_arr[i]);
        free(header_arr[i]);
    }

    fclose(output_file);
    fclose(index_file);
    fclose(query_file);

    free(query_arr);
    free(header_arr);
    free(text_buffer);
    free(genome);
    free(suffix_arr);

    free(bwt);
    free(tally_a);
    free(tally_c);
    free(tally_g);
    free(tally_t);


    printf("Done running querysa\n^-^\n");

    return 0;
}