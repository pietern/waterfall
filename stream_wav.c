#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <arpa/inet.h>
#include <alloca.h>
#include <string.h>

#include "stream_raw.h"

typedef struct _wav_header_s _wav_header_t;

struct _wav_header_s {
    char riff[4];
    uint32_t file_size;
    char wave[4];

    char fmt[4];
    uint32_t fmt_len;

    uint16_t type;
    uint16_t channels;
    uint32_t sample_rate;
    uint32_t _ignore_1;
    uint16_t _ignore_2;
    uint16_t bits_per_sample;

    char data[4];
    uint32_t data_len;

    // Following fields not included in file header
    uint16_t bytes_per_sample;
};

static int8_t
stream_wav__init(__attribute__((unused)) const stream_handler_t *h,
                 stream_t *s)
{
    _wav_header_t *header;
    int8_t rv;

    header = malloc(sizeof(*header));
    rv = read(s->fd, header, 44);
    if (rv < 0) {
        return rv;
    }

    if (strncmp(header->riff, "RIFF", 4) != 0) {
        return -1;
    }

    if (strncmp(header->wave, "WAVE", 4) != 0) {
        return -1;
    }

    if (header->channels != 2) {
        return -1;
    }

    header->bytes_per_sample = (2 * (header->bits_per_sample / 8));

    s->sample_rate = header->sample_rate;
    s->n_samples = (s->stat.st_size - sizeof(*header)) / header->bytes_per_sample;
    s->data = header;

    return 0;
}

static uint64_t
stream_wav__to_byte_offset(__attribute__((unused)) const stream_handler_t *h,
                           stream_t *s,
                           uint64_t offset)
{
    const _wav_header_t *header;
    header = (const _wav_header_t *) s->data;
    return sizeof(_wav_header_t) + 2 * (header->bits_per_sample / 8) * offset;
}

static uint64_t
stream_wav__to_byte_length(__attribute__((unused)) const stream_handler_t *h,
                           stream_t *s,
                           uint64_t length)
{
    const _wav_header_t *header;
    header = (const _wav_header_t *) s->data;
    return 2 * (header->bits_per_sample / 8) * length;
}

static int8_t
stream_wav__read(const stream_handler_t *h,
                 stream_t *s,
                 uint64_t offset,
                 uint64_t length,
                 fftw_complex *out)
{
    uint64_t byte_offset = stream_wav__to_byte_offset(h, s, offset);
    uint64_t byte_length = stream_wav__to_byte_length(h, s, length);
    off_t rv;
    const _wav_header_t *header;
    uint64_t i;

    rv = lseek(s->fd, byte_offset, SEEK_SET);
    if (rv < 0) {
        return rv;
    }

    header = (const _wav_header_t *) s->data;
    if (header->bits_per_sample == 16) {
        int16_t *in;

        in = alloca(byte_length);
        rv = read(s->fd, in, byte_length);
        if (rv < 0) {
            return rv;
        }

        for (i = 0; i < length; i++) {
            out[i][0] = (-(float) in[i * 2 + 0]) / INT16_MIN;
            out[i][1] = (-(float) in[i * 2 + 1]) / INT16_MIN;
        }
    } else {
        assert(0 && "invalid sample size");
    }

    return 0;
}

stream_handler_t stream_wav_handler = {
    .init = &stream_wav__init,
    .to_byte_offset = &stream_wav__to_byte_offset,
    .to_byte_length = &stream_wav__to_byte_length,
    .read = &stream_wav__read,
};
