#include <stdlib.h>
#include <sys/mman.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include <alloca.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/fcntl.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>

#include "stream.h"

stream_t *stream_open(const char *file, const stream_handler_t *handler) {
    stream_t *s;
    int rv;

    s = calloc(1, sizeof(*s));
    s->fd = -1;
    s->file = file;
    s->handler = handler;
    s->sample_rate = -1;
    s->center_freq = -1;

    s->fd = open(file, O_RDONLY);
    if (s->fd < 0) {
        goto error;
    }

    rv = fstat(s->fd, &s->stat);
    if (rv < 0) {
        goto error;
    }

#if _XOPEN_SOURCE >= 600 || _POSIX_C_SOURCE >= 200112L
    rv = posix_fadvise(s->fd, 0, s->stat.st_size, POSIX_FADV_RANDOM);
    if (rv < 0) {
        goto error;
    }
#endif

    rv = s->handler->init(s->handler, s);
    if (rv < 0) {
        goto error;
    }

    return s;

 error:
    stream_close(s);
    return NULL;
}

void stream_close(stream_t *s) {
    if (s != NULL) {
        if (s->fd >= 0) {
            close(s->fd);
        }
        free(s);
    }
}

void stream_fftw_init(stream_t *s) {
    int i;

    for (i = 0; i < 20; i++) {
        s->fftw_in[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (1 << i));
        s->fftw_out[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (1 << i));
        s->fftw_plan[i] = fftw_plan_dft_1d(1 << i, s->fftw_in[i], s->fftw_out[i], FFTW_FORWARD, FFTW_ESTIMATE);
        s->fftw_n[i] = 1 << i;
    }
}

void stream_fftw_free(stream_t *s) {
    int i;

    for (i = 0; i < 20; i++) {
        fftw_destroy_plan(s->fftw_plan[i]);
        fftw_free(s->fftw_in[i]);
        fftw_free(s->fftw_out[i]);
    }
}

static void stream_fft_process(stream_t *s, uint8_t scale, double *dst) {
    uint16_t i;
    uint16_t h = s->fftw_n[scale] / 2;
    uint16_t p;
    double mag;

    for (i = 0; i < s->fftw_n[scale]; i++) {
        // Shift FFT
        if (i < h) {
            p = h + i;
        } else {
            p = i - h;
        }

        // Compute and normalize magnitude
        mag = (s->fftw_out[scale][p][0] * s->fftw_out[scale][p][0] +
               s->fftw_out[scale][p][1] * s->fftw_out[scale][p][1]);
        mag = mag / (s->fftw_n[scale] * s->fftw_n[scale]);

        // Convert to dB
        mag = 10.0f * log10f(mag);

        // Accumulate magnitude in destination buffer
        dst[i] += mag;
    }
}

// From: https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
static uint8_t log2i(uint32_t v) {
    static const unsigned int b[] = {
        0xAAAAAAAA,
        0xCCCCCCCC,
        0xF0F0F0F0,
        0xFF00FF00,
        0xFFFF0000,
    };
    uint8_t i;
    uint8_t r;

    r = (v & b[0]) != 0;
    for (i = 4; i > 0; i--) // unroll for speed...
    {
        r |= ((v & b[i]) != 0) << i;
    }

    return r;
}

void stream_fft(stream_t *s, uint64_t t, double *d, uint32_t len) {
    int rv;
    int scale;

    if ((len & (len - 1)) != 0x0) {
        assert(0 && "not a power of 2");
    }

    scale = log2i(len);
    rv = s->handler->read(s->handler, s, t, s->fftw_n[scale], s->fftw_in[scale]);
    if (rv < 0) {
        perror("read");
        assert(0);
    }

    fftw_execute(s->fftw_plan[scale]);
    stream_fft_process(s, scale, d);
}


static void stream__fadvise(stream_t *s, int8_t advise, uint64_t start, uint64_t len) {
    uint64_t offset;
    uint64_t bytes;
    int8_t rv;

    offset = s->handler->to_byte_offset(s->handler, s, start);
    bytes = s->handler->to_byte_length(s->handler, s, len);
    rv = posix_fadvise(s->fd, offset, bytes, advise);
    if (rv < 0) {
        perror("posix_fadvise");
        assert(0);
    }
}

void stream_fadvise_willneed(stream_t *s, uint64_t start, uint64_t len) {
    stream__fadvise(s, POSIX_FADV_WILLNEED, start, len);
}

void stream_fadvise_dontneed(stream_t *s, uint64_t start, uint64_t len) {
    stream__fadvise(s, POSIX_FADV_DONTNEED, start, len);
}
