#ifndef _STREAM_H
#define _STREAM_H

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdint.h>
#include <fftw3.h>

typedef struct stream_s stream_t;

typedef struct stream_handler_s stream_handler_t;

struct stream_handler_s {
    int8_t (*init)(const stream_handler_t *, stream_t *);
    uint64_t (*to_byte_offset)(const stream_handler_t *, stream_t *, uint64_t);
    uint64_t (*to_byte_length)(const stream_handler_t *, stream_t *, uint64_t);
    int8_t (*read)(const stream_handler_t *, stream_t *, uint64_t, uint64_t, fftw_complex *);

    // Private data for handler.
    // Keep in mind that a single handler can be used for multiple
    // streams, so any data specific to the stream should be stored
    // in the data field of the stream itself.
    void *data;
};

struct stream_s {
    const char *file;
    const stream_handler_t *handler;

    int fd;
    struct stat stat;
    void *data;

    int32_t sample_rate;
    int32_t center_freq;
    int64_t n_samples;

    fftw_complex *fftw_in[20];
    fftw_complex *fftw_out[20];
    fftw_plan fftw_plan[20];
    uint16_t fftw_n[20];
};

stream_t *stream_open(const char *file, const stream_handler_t *handler);

void stream_close(stream_t *s);

void stream_fftw_init(stream_t *s);

void stream_fftw_free(stream_t *s);

void stream_fft(stream_t *s, uint64_t t, double *d, uint32_t len);

void stream_fadvise_willneed(stream_t *s, uint64_t start, uint64_t len);

void stream_fadvise_dontneed(stream_t *s, uint64_t start, uint64_t len);

#endif
