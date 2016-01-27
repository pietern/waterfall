#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/mman.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <strings.h>

#include "stream.h"
#include "window.h"
#include "wedge.h"
#include "color.h"

#include "stream_raw.h"
#include "stream_wav.h"

static struct {
    char *type;
    int32_t sample_rate;
    int32_t center_freq;
} *global;

void usage(__attribute__((unused)) int argc, char **argv) {
    fprintf(stderr, "Usage: %s [OPTION]... FILE\n", argv[0]);
    fprintf(stderr, "Analyze frequency spectrum of signal in FILE\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -t TYPE    File type of FILE [wav,s8,s16,s32,f32]\n");
    fprintf(stderr, "  -r RATE    Sample rate of signal\n");
    fprintf(stderr, "  -c FREQ    Center frequency of signal\n");
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
}

static void parse(int argc, char **argv) {
    int opt;

    if (argc < 2) {
        usage(argc, argv);
    }

    if (strcmp(argv[1], "--help") == 0) {
        usage(argc, argv);
    }

    global = calloc(sizeof(*global), 1);
    global->type = NULL;
    global->sample_rate = -1;
    global->center_freq = -1;

    while ((opt = getopt(argc, argv, "ht:r:c:")) != -1) {
        switch (opt) {
        case 'h':
            usage(argc, argv);
            break;
        case 't':
            global->type = strdup(optarg);
            break;
        case 'r':
            global->sample_rate = atoi(optarg);
            break;
        case 'c':
            global->center_freq = atoi(optarg);
            break;
        default:
            fprintf(stderr, "Try '%s --help' for more information.\n", argv[0]);
            exit(EXIT_FAILURE);
        }
    }
}

void detect_handler(const char *file) {
    const char *ext;

    ext = strrchr(file, '.');
    if (ext == NULL) {
        return;
    }

    if (strcasecmp(ext + 1, "wav") == 0) {
        global->type = "wav";
    }
}

int main(int argc, char **argv) {
    const char *file;
    stream_handler_t *handler;
    stream_t *stream;

    parse(argc, argv);

    if (optind < argc) {
        file = argv[optind];
    } else {
        fprintf(stderr, "%s: please specify file\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    // If type is not specified as argument,
    // try to derive it from the file name
    if (global->type == NULL) {
        detect_handler(file);
    }

    // Initialize handler from file type
    if (global->type != NULL) {
        if (strcasecmp(global->type, "wav") == 0) {
            handler = &stream_wav_handler;
        } else if (strcasecmp(global->type, "s8") == 0) {
            handler = &stream_s8_handler;
        } else if (strcasecmp(global->type, "s16") == 0) {
            handler = &stream_s16_handler;
        } else if (strcasecmp(global->type, "s32") == 0) {
            handler = &stream_s32_handler;
        } else if (strcasecmp(global->type, "f32") == 0) {
            handler = &stream_f32_handler;
        } else {
            fprintf(stderr, "%s: invalid file type '%s'\n", argv[0], global->type);
            exit(EXIT_FAILURE);
        }
    } else {
        fprintf(stderr, "%s: please specify file type\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    stream = stream_open(file, handler);
    if (stream == NULL) {
        fprintf(stderr, "%s: error opening stream\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    // Initialize sample rate if not yet initialized by the handler
    if (stream->sample_rate <= 0) {
        if (global->sample_rate <= 0) {
            fprintf(stderr, "%s: please specify sample rate\n", argv[0]);
            exit(EXIT_FAILURE);
        }

        stream->sample_rate = global->sample_rate;
    }

    // Initialize center frequency
    if (stream->center_freq <= 0) {
        if (global->center_freq <= 0) {
            stream->center_freq = 0;
        } else {
            stream->center_freq = global->center_freq;
        }
    }

    stream_fftw_init(stream);

    wedge_processing_init();

    color = &colormap_rainbow[0];

    window_run(stream);

    stream_close(stream);

    exit(EXIT_SUCCESS);
}
