#include <sys/types.h>
#include <unistd.h>
#include <assert.h>
#include <stdlib.h>

#include "stream_raw.h"

typedef struct _raw_s _raw_t;

struct _raw_s {
    // Bytes per sample component (in-phase or quadrature)
    uint8_t bytes;

    // Function to process arbitrary data into complex samples
    void (*process)(const void *, fftw_complex *, uint64_t);
};

static int8_t
stream_raw__init(const stream_handler_t *h,
                 stream_t *s)
{
    _raw_t *r = (_raw_t *) h->data;
    s->n_samples = s->stat.st_size / (2 * r->bytes);
    return 0;
}

static uint64_t
stream_raw__to_byte_offset(const stream_handler_t *h,
                           __attribute__((unused)) stream_t *s,
                           uint64_t offset)
{
    _raw_t *r = (_raw_t *) h->data;
    return 2 * r->bytes * offset;
}

static uint64_t
stream_raw__to_byte_length(const stream_handler_t *h,
                           __attribute__((unused)) stream_t *s,
                           uint64_t length)
{
    const _raw_t *r = (const _raw_t *) h->data;
    return 2 * r->bytes * length;
}

static int8_t
stream_raw__read(const stream_handler_t *h,
                 stream_t *s,
                 uint64_t offset,
                 uint64_t length,
                 fftw_complex *out)
{
    const _raw_t *r = (const _raw_t *) h->data;
    uint64_t byte_offset = stream_raw__to_byte_offset(h, s, offset);
    uint64_t byte_length = stream_raw__to_byte_length(h, s, length);
    off_t rv;
    void *in;

    rv = lseek(s->fd, byte_offset, SEEK_SET);
    if (rv < 0) {
        return rv;
    }

    in = malloc(byte_length);
    rv = read(s->fd, in, byte_length);
    if (rv < 0) {
        return rv;
    }

    r->process(in, out, length);
    free(in);
    return 0;
}

static void
_raw_s8_process(const void *data,
                fftw_complex *out,
                uint64_t length)
{
    const int8_t *in = (const int8_t *) data;
    uint64_t i;

    for (i = 0; i < length; i++) {
        out[i][0] = (-(float) in[i * 2 + 0]) / INT8_MIN;
        out[i][1] = (-(float) in[i * 2 + 1]) / INT8_MIN;
    }
}

static _raw_t raw_s8 = {
    .bytes = 1,
    .process = &_raw_s8_process,
};

stream_handler_t stream_s8_handler = {
    .init = &stream_raw__init,
    .to_byte_offset = &stream_raw__to_byte_offset,
    .to_byte_length = &stream_raw__to_byte_length,
    .read = &stream_raw__read,
    .data = &raw_s8,
};

static void
_raw_s16_process(const void *data,
                 fftw_complex *out,
                 uint64_t length)
{
    const int16_t *in = (const int16_t *) data;
    uint64_t i;

    for (i = 0; i < length; i++) {
        out[i][0] = (-(float) in[i * 2 + 0]) / INT16_MIN;
        out[i][1] = (-(float) in[i * 2 + 1]) / INT16_MIN;
    }
}

static _raw_t raw_s16 = {
    .bytes = 2,
    .process = &_raw_s16_process,
};

stream_handler_t stream_s16_handler = {
    .init = &stream_raw__init,
    .to_byte_offset = &stream_raw__to_byte_offset,
    .to_byte_length = &stream_raw__to_byte_length,
    .read = &stream_raw__read,
    .data = &raw_s16,
};

static void
_raw_s32_process(const void *data,
                 fftw_complex *out,
                 uint64_t length)
{
    const int32_t *in = (const int32_t *) data;
    uint64_t i;

    for (i = 0; i < length; i++) {
        out[i][0] = (-(float) in[i * 2 + 0]) / INT32_MIN;
        out[i][1] = (-(float) in[i * 2 + 1]) / INT32_MIN;
    }
}

static _raw_t raw_s32 = {
    .bytes = 4,
    .process = &_raw_s32_process,
};

stream_handler_t stream_s32_handler = {
    .init = &stream_raw__init,
    .to_byte_offset = &stream_raw__to_byte_offset,
    .to_byte_length = &stream_raw__to_byte_length,
    .read = &stream_raw__read,
    .data = &raw_s32,
};

static void
_raw_f32_process(const void *data,
                 fftw_complex *out,
                 uint64_t length)
{
    const float *in = (const float *) data;
    uint64_t i;

    for (i = 0; i < length; i++) {
        out[i][0] = in[i * 2 + 0];
        out[i][1] = in[i * 2 + 1];
    }
}

static _raw_t raw_f32 = {
    .bytes = 4,
    .process = &_raw_f32_process,
};

stream_handler_t stream_f32_handler = {
    .init = &stream_raw__init,
    .to_byte_offset = &stream_raw__to_byte_offset,
    .to_byte_length = &stream_raw__to_byte_length,
    .read = &stream_raw__read,
    .data = &raw_f32,
};
