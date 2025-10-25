/*
  huffman_serial.c
  Serial Huffman compressor/decompressor (single-file).
  Usage:
    ./huffman_serial c input.bin out.huff
    ./huffman_serial d out.huff recovered.bin
*/

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <errno.h>

#define SYMBOLS 256
#define MAGIC 0x48465546U /* "HFUF" */

typedef struct Node {
    uint64_t freq;
    int symbol; /* -1 for internal nodes */
    struct Node *left, *right;
} Node;

typedef struct {
    uint8_t *bits;   /* bits stored MSB-first in bytes */
    size_t bitlen;   /* number of bits */
} Code;

/* ------------------ utilities ------------------ */

static void *xmalloc(size_t n) {
    void *p = malloc(n);
    if (!p) { perror("malloc"); exit(1); }
    return p;
}

static uint64_t now_ns() {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (uint64_t)t.tv_sec * 1000000000ULL + t.tv_nsec;
}

/* ------------------ min-heap for nodes ------------------ */

typedef struct {
    Node **a;
    int n;
    int cap;
} MinHeap;

MinHeap *heap_create(int cap) {
    MinHeap *h = xmalloc(sizeof(MinHeap));
    h->a = xmalloc(sizeof(Node*) * (cap+4));
    h->n = 0; h->cap = cap+4;
    return h;
}
void heap_push(MinHeap *h, Node *node) {
    int i = h->n++;
    h->a[i] = node;
    while (i > 0) {
        int p = (i - 1) >> 1;
        if (h->a[p]->freq <= h->a[i]->freq) break;
        Node *tmp = h->a[p]; h->a[p] = h->a[i]; h->a[i] = tmp;
        i = p;
    }
}
Node* heap_pop(MinHeap *h) {
    if (h->n == 0) return NULL;
    Node *ret = h->a[0];
    h->a[0] = h->a[--h->n];
    int i = 0;
    for (;;) {
        int l = i*2 + 1, r = i*2 + 2, s = i;
        if (l < h->n && h->a[l]->freq < h->a[s]->freq) s = l;
        if (r < h->n && h->a[r]->freq < h->a[s]->freq) s = r;
        if (s == i) break;
        Node *tmp = h->a[i]; h->a[i] = h->a[s]; h->a[s] = tmp;
        i = s;
    }
    return ret;
}
void heap_free(MinHeap *h) { free(h->a); free(h); }

/* ------------------ Huffman tree & code generation ------------------ */

Node *node_new(uint64_t f, int sym, Node *l, Node *r) {
    Node *n = xmalloc(sizeof(Node));
    n->freq = f; n->symbol = sym; n->left = l; n->right = r;
    return n;
}

Node* build_tree(uint64_t freq[SYMBOLS]) {
    MinHeap *h = heap_create(SYMBOLS);
    for (int i = 0; i < SYMBOLS; ++i) {
        if (freq[i] > 0) heap_push(h, node_new(freq[i], i, NULL, NULL));
    }
    if (h->n == 0) {
        /* empty input: create single node for symbol 0 */
        heap_push(h, node_new(1, 0, NULL, NULL));
    }
    while (h->n > 1) {
        Node *a = heap_pop(h);
        Node *b = heap_pop(h);
        Node *c = node_new(a->freq + b->freq, -1, a, b);
        heap_push(h, c);
    }
    Node *root = heap_pop(h);
    heap_free(h);
    return root;
}

void free_tree(Node *n) {
    if (!n) return;
    free_tree(n->left);
    free_tree(n->right);
    free(n);
}

/* generate codes by traversing tree; store bits in dynamic arrays */
void gen_codes_rec(Node *n, Code codes[SYMBOLS], uint8_t *path, int depth) {
    if (!n) return;
    if (n->symbol >= 0) {
        Code *c = &codes[n->symbol];
        c->bitlen = depth;
        size_t bytes = (depth + 7) >> 3;
        c->bits = xmalloc(bytes ? bytes : 1);
        memset(c->bits, 0, bytes ? bytes : 1);
        for (int i = 0; i < depth; ++i) {
            if (path[i]) c->bits[i>>3] |= (1 << (7 - (i & 7)));
        }
        return;
    }
    path[depth] = 0;
    gen_codes_rec(n->left, codes, path, depth+1);
    path[depth] = 1;
    gen_codes_rec(n->right, codes, path, depth+1);
}

void gen_codes(Node *root, Code codes[SYMBOLS]) {
    for (int i = 0; i < SYMBOLS; ++i) { codes[i].bits = NULL; codes[i].bitlen = 0; }
    uint8_t *path = xmalloc(512);
    gen_codes_rec(root, codes, path, 0);
    free(path);
}

void free_codes(Code codes[SYMBOLS]) {
    for (int i = 0; i < SYMBOLS; ++i) if (codes[i].bits) free(codes[i].bits);
}

/* ------------------ bit writer ------------------ */

typedef struct {
    uint8_t *buf;
    size_t cap;
    size_t size; /* bytes fully used */
    uint8_t bitpos; /* 0..7 next free bit index in current byte */
} BitWriter;

BitWriter* bw_create(size_t initial) {
    BitWriter *w = xmalloc(sizeof(BitWriter));
    w->cap = initial ? initial : 1024;
    w->buf = xmalloc(w->cap);
    w->size = 0; w->bitpos = 0;
    memset(w->buf, 0, w->cap);
    return w;
}
void bw_ensure(BitWriter *w, size_t more_bits) {
    size_t need_bytes = ((w->bitpos + more_bits) + 7) >> 3;
    if (w->size + need_bytes >= w->cap) {
        size_t newcap = (w->cap * 3) / 2 + need_bytes + 16;
        uint8_t *nb = realloc(w->buf, newcap);
        if (!nb) { perror("realloc"); exit(1); }
        memset(nb + w->cap, 0, newcap - w->cap);
        w->buf = nb; w->cap = newcap;
    }
}
void bw_put_bits(BitWriter *w, const uint8_t *bits, size_t bitlen) {
    bw_ensure(w, bitlen);
    for (size_t i = 0; i < bitlen; ++i) {
        int val = (bits[i>>3] >> (7 - (i & 7))) & 1;
        if (w->bitpos == 0 && w->size == 0) {
            /* nothing yet, make sure first byte exists (we keep buf zeroed) */
        }
        if (w->bitpos == 0 && w->size > 0 && ((w->size) >= w->cap)) {
            /* expanded earlier */
        }
        if (w->bitpos == 0 && w->size >= w->cap) { /* unlikely */ bw_ensure(w, 8); }
        if (w->bitpos == 0) {
            /* start new byte */
            w->buf[w->size] = 0;
        }
        w->buf[w->size] |= (val << (7 - w->bitpos));
        w->bitpos++;
        if (w->bitpos == 8) { w->bitpos = 0; w->size++; }
    }
}
void bw_finish(BitWriter *w) {
    if (w->bitpos != 0) w->size++;
}
void bw_free(BitWriter *w) {
    if (!w) return;
    free(w->buf);
    free(w);
}

/* ------------------ bit reader helpers ------------------ */

typedef struct {
    const uint8_t *buf;
    size_t byte_len;
    uint64_t bit_len;
} BitStream;

/* ------------------ encode/decode ------------------ */

/* encode input buffer into a bit buffer using codes[] */
void encode_to_bitstream(const uint8_t *inbuf, size_t inlen, Code codes[SYMBOLS],
                         uint8_t **out_bytes, size_t *out_byte_len, uint64_t *out_bitlen) {
    BitWriter *bw = bw_create(inlen / 2 + 16);
    uint64_t total_bits = 0;
    for (size_t i = 0; i < inlen; ++i) {
        uint8_t s = inbuf[i];
        bw_put_bits(bw, codes[s].bits, codes[s].bitlen);
        total_bits += codes[s].bitlen;
    }
    bw_finish(bw);
    *out_bytes = bw->buf;
    *out_byte_len = bw->size;
    *out_bitlen = total_bits;
    free(bw); /* NOTE: freed but we 'stole' bw->buf; therefore this free is incorrect - fix below */
}

/* We need a version that returns the buffer without freeing it accidentally; rewrite properly */
void encode_to_bitstream_safe(const uint8_t *inbuf, size_t inlen, Code codes[SYMBOLS],
                         uint8_t **out_bytes, size_t *out_byte_len, uint64_t *out_bitlen) {
    BitWriter *bw = bw_create(inlen / 2 + 16);
    uint64_t total_bits = 0;
    for (size_t i = 0; i < inlen; ++i) {
        uint8_t s = inbuf[i];
        bw_put_bits(bw, codes[s].bits, codes[s].bitlen);
        total_bits += codes[s].bitlen;
    }
    bw_finish(bw);
    *out_bytes = bw->buf;
    *out_byte_len = bw->size;
    *out_bitlen = total_bits;
    /* don't free bw->buf here; caller will free *out_bytes */
    free(bw);
}

/* decode bitstream using Huffman tree; stop after producing expected_out bytes (original size) */
uint8_t* decode_from_bitstream(const uint8_t *bits, uint64_t bitlen, Node *root, size_t expected_out, size_t *produced_out) {
    uint8_t *out = xmalloc(expected_out ? expected_out : 1);
    size_t outpos = 0;
    Node *cur = root;
    for (uint64_t i = 0; i < bitlen; ++i) {
        int val = (bits[i>>3] >> (7 - (i & 7))) & 1;
        cur = val ? cur->right : cur->left;
        if (!cur) { fprintf(stderr, "Decoding error: traversed to NULL\n"); free(out); *produced_out = 0; return NULL; }
        if (cur->left == NULL && cur->right == NULL) {
            /* leaf */
            if (outpos >= expected_out) {
                /* extra safety: should not happen */
                break;
            }
            out[outpos++] = (uint8_t)cur->symbol;
            cur = root;
        }
    }
    *produced_out = outpos;
    return out;
}

/* ------------------ file IO helpers ------------------ */

static void write_u32(FILE *f, uint32_t v) { if (fwrite(&v, sizeof(v), 1, f) != 1) perror("write_u32"); }
static void write_u64(FILE *f, uint64_t v) { if (fwrite(&v, sizeof(v), 1, f) != 1) perror("write_u64"); }
static uint32_t read_u32(FILE *f) { uint32_t v=0; if (fread(&v, sizeof(v), 1, f) != 1) { if (!feof(f)) perror("read_u32"); } return v; }
static uint64_t read_u64(FILE *f) { uint64_t v=0; if (fread(&v, sizeof(v), 1, f) != 1) { if (!feof(f)) perror("read_u64"); } return v; }

uint8_t* read_full(const char *path, size_t *out_size) {
    FILE *f = fopen(path, "rb");
    if (!f) { perror("fopen"); return NULL; }
    struct stat st;
    if (stat(path, &st) == 0) *out_size = st.st_size;
    else { fseek(f, 0, SEEK_END); *out_size = ftell(f); fseek(f, 0, SEEK_SET); }
    uint8_t *buf = xmalloc(*out_size ? *out_size : 1);
    if (*out_size) {
        size_t r = fread(buf, 1, *out_size, f);
        if (r != *out_size) { fprintf(stderr, "Short read\n"); }
    }
    fclose(f);
    return buf;
}

/* ------------------ compress / decompress ------------------ */

int compress_file(const char *inpath, const char *outpath) {
    size_t inlen;
    uint8_t *in = read_full(inpath, &inlen);
    if (!in) return 1;

    uint64_t freq[SYMBOLS]; memset(freq,0,sizeof(freq));
    for (size_t i = 0; i < inlen; ++i) freq[in[i]]++;

    Node *root = build_tree(freq);
    Code codes[SYMBOLS];
    gen_codes(root, codes);

    /* encode entire file */
    uint8_t *bitbytes = NULL;
    size_t bitbytes_len = 0;
    uint64_t bitlen = 0;
    encode_to_bitstream_safe(in, inlen, codes, &bitbytes, &bitbytes_len, &bitlen);

    /* write output */
    FILE *fo = fopen(outpath, "wb");
    if (!fo) { perror("fopen out"); free(bitbytes); free(in); free_codes(codes); free_tree(root); return 1; }
    write_u32(fo, MAGIC);
    write_u64(fo, (uint64_t)inlen);
    for (int i = 0; i < SYMBOLS; ++i) write_u64(fo, freq[i]);
    write_u64(fo, bitlen);
    write_u64(fo, (uint64_t)bitbytes_len);
    if (bitbytes_len) {
        size_t w = fwrite(bitbytes, 1, bitbytes_len, fo);
        if (w != bitbytes_len) { perror("fwrite bitbytes"); fclose(fo); free(bitbytes); free(in); free_codes(codes); free_tree(root); return 1; }
    }
    fclose(fo);

    free(bitbytes);
    free_codes(codes);
    free_tree(root);
    free(in);

    return 0;
}

int decompress_file(const char *inpath, const char *outpath) {
    FILE *fi = fopen(inpath, "rb");
    if (!fi) { perror("fopen in"); return 1; }
    uint32_t magic = read_u32(fi);
    if (magic != MAGIC) { fprintf(stderr, "Not a supported file (bad magic)\n"); fclose(fi); return 1; }
    uint64_t orig_size = read_u64(fi);
    uint64_t freq[SYMBOLS];
    for (int i = 0; i < SYMBOLS; ++i) freq[i] = read_u64(fi);
    uint64_t bitlen = read_u64(fi);
    uint64_t bitbytes_len = read_u64(fi);
    uint8_t *bitbytes = NULL;
    if (bitbytes_len) {
        bitbytes = xmalloc(bitbytes_len);
        size_t r = fread(bitbytes, 1, bitbytes_len, fi);
        if (r != bitbytes_len) { fprintf(stderr, "Short read of bitstream\n"); free(bitbytes); fclose(fi); return 1; }
    }
    fclose(fi);

    Node *root = build_tree(freq);
    size_t produced = 0;
    uint8_t *out = NULL;
    if (bitlen == 0 || bitbytes_len == 0) {
        /* nothing encoded -> produce empty output */
        out = xmalloc(1);
        produced = 0;
    } else {
        out = decode_from_bitstream(bitbytes, bitlen, root, (size_t)orig_size, &produced);
    }

    /* write out */
    FILE *fo = fopen(outpath, "wb");
    if (!fo) { perror("fopen out"); free(bitbytes); free(out); free_tree(root); return 1; }
    if (produced) {
        size_t w = fwrite(out, 1, produced, fo);
        if (w != produced) { perror("fwrite out"); fclose(fo); free(bitbytes); free(out); free_tree(root); return 1; }
    }
    fclose(fo);

    free(bitbytes);
    free(out);
    free_tree(root);

    return 0;
}

/* ------------------ main ------------------ */

void print_usage(const char *p) {
    fprintf(stderr, "Usage:\n  %s c input output    # compress\n  %s d input output    # decompress\n", p, p);
}

int main(int argc, char **argv) {
    if (argc != 4) { print_usage(argv[0]); return 1; }
    char op = argv[1][0];
    const char *inpath = argv[2];
    const char *outpath = argv[3];

    if (op == 'c') {
        uint64_t t0 = now_ns();
        int r = compress_file(inpath, outpath);
        uint64_t t1 = now_ns();
        if (r == 0) printf("Compression done in %.6f s\n", (t1 - t0) / 1e9);
        return r;
    } else if (op == 'd') {
        uint64_t t0 = now_ns();
        int r = decompress_file(inpath, outpath);
        uint64_t t1 = now_ns();
        if (r == 0) printf("Decompression done in %.6f s\n", (t1 - t0) / 1e9);
        return r;
    } else {
        print_usage(argv[0]);
        return 1;
    }
}
