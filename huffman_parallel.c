#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <pthread.h>
#include <sys/stat.h>

#define DEFAULT_CHUNK_SIZE (1<<20) // 1MB
#define DEFAULT_THREADS 4

typedef struct {
    unsigned char *data;
    size_t size;
} Chunk;

// Placeholder for your Huffman compress function
void huffman_compress_chunk(unsigned char *input, size_t size, unsigned char **out, size_t *out_size) {
    // TODO: replace with actual Huffman compression
    *out = malloc(size); 
    memcpy(*out, input, size);
    *out_size = size;
}

// Placeholder for your Huffman decompress function
void huffman_decompress_chunk(unsigned char *input, size_t size, unsigned char **out, size_t *out_size) {
    // TODO: replace with actual Huffman decompression
    *out = malloc(size); 
    memcpy(*out, input, size);
    *out_size = size;
}

int main(int argc, char *argv[]) {
    if(argc < 5) {
        printf("Usage: %s c|d serial|openmp|pthreads input output [--chunk N] [--threads T]\n", argv[0]);
        return 1;
    }

    char mode = argv[1][0]; // 'c' or 'd'
    char *parallel_mode = argv[2]; // serial|openmp|pthreads
    char *infile = argv[3];
    char *outfile = argv[4];

    size_t chunk_size = DEFAULT_CHUNK_SIZE;
    int nthreads = DEFAULT_THREADS;

    // Parse optional args
    for(int i=5;i<argc;i++) {
        if(strcmp(argv[i], "--chunk")==0 && i+1<argc) chunk_size = atol(argv[++i]);
        if(strcmp(argv[i], "--threads")==0 && i+1<argc) nthreads = atoi(argv[++i]);
    }

    // Read input file
    FILE *fin = fopen(infile,"rb");
    if(!fin) { perror("fopen"); return 1; }

    fseek(fin, 0, SEEK_END);
    size_t fsize = ftell(fin);
    fseek(fin, 0, SEEK_SET);

    unsigned char *buffer = malloc(fsize);
    fread(buffer, 1, fsize, fin);
    fclose(fin);

    size_t n_chunks = (fsize + chunk_size - 1)/chunk_size;
    Chunk *chunks = malloc(sizeof(Chunk)*n_chunks);

    for(size_t i=0;i<n_chunks;i++) {
        size_t sz = (i==n_chunks-1)? (fsize - i*chunk_size) : chunk_size;
        chunks[i].data = buffer + i*chunk_size;
        chunks[i].size = sz;
    }

    unsigned char **out_chunks = malloc(sizeof(unsigned char*)*n_chunks);
    size_t *out_sizes = malloc(sizeof(size_t)*n_chunks);

    double t0 = omp_get_wtime();

    if(strcmp(parallel_mode,"serial")==0) {
        for(size_t i=0;i<n_chunks;i++) {
            if(mode=='c') huffman_compress_chunk(chunks[i].data, chunks[i].size, &out_chunks[i], &out_sizes[i]);
            else huffman_decompress_chunk(chunks[i].data, chunks[i].size, &out_chunks[i], &out_sizes[i]);
        }
    } else if(strcmp(parallel_mode,"openmp")==0) {
        omp_set_num_threads(nthreads);
        #pragma omp parallel for
        for(size_t i=0;i<n_chunks;i++) {
            if(mode=='c') huffman_compress_chunk(chunks[i].data, chunks[i].size, &out_chunks[i], &out_sizes[i]);
            else huffman_decompress_chunk(chunks[i].data, chunks[i].size, &out_chunks[i], &out_sizes[i]);
        }
    } else if(strcmp(parallel_mode,"pthreads")==0) {
        // Simple pthread wrapper
        pthread_t *threads = malloc(sizeof(pthread_t)*nthreads);
        typedef struct { size_t start, end; } ThreadData;
        ThreadData *td = malloc(sizeof(ThreadData)*nthreads);

        size_t per_thread = (n_chunks + nthreads - 1)/nthreads;
        void *thread_func(void *arg) {
            ThreadData *d = (ThreadData*)arg;
            for(size_t i=d->start;i<d->end && i<n_chunks;i++) {
                if(mode=='c') huffman_compress_chunk(chunks[i].data, chunks[i].size, &out_chunks[i], &out_sizes[i]);
                else huffman_decompress_chunk(chunks[i].data, chunks[i].size, &out_chunks[i], &out_sizes[i]);
            }
            return NULL;
        }

        for(int i=0;i<nthreads;i++) {
            td[i].start = i*per_thread;
            td[i].end = td[i].start + per_thread;
            pthread_create(&threads[i], NULL, thread_func, &td[i]);
        }
        for(int i=0;i<nthreads;i++) pthread_join(threads[i], NULL);
        free(threads); free(td);
    } else {
        printf("Unknown mode\n"); return 1;
    }

    double t1 = omp_get_wtime();
    printf("%sion (%s) took %f s\n", (mode=='c')?"Compress":"Decompress", parallel_mode, t1-t0);

    // Write output file
    FILE *fout = fopen(outfile,"wb");
    if(!fout) { perror("fopen"); return 1; }
    for(size_t i=0;i<n_chunks;i++) {
        fwrite(out_chunks[i],1,out_sizes[i],fout);
        free(out_chunks[i]);
    }
    fclose(fout);

    free(out_chunks); free(out_sizes); free(chunks); free(buffer);
    return 0;
}
