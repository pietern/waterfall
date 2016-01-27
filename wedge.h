#ifndef _WEDGE_H
#define _WEDGE_H

#include <stdint.h>
#include <gtk/gtk.h>

#include "stream.h"

#define WEDGE_QUEUED 1
#define WEDGE_PROCESSING 2
#define WEDGE_DONE 3

typedef struct wedge_s wedge_t;

struct wedge_s {
    volatile int64_t status;

    // Callback number when this wedge was most recently used
    uint64_t draw_cb;

    stream_t *stream;

    int width;
    int height;
    int stride;

    int sx;
    int sy;
    int64_t start;
    int64_t end;
    int64_t skip;

    unsigned char *data;
    cairo_surface_t *surface;

    void (*callback_fn)(wedge_t *, void *);
    void *callback_data;
};

wedge_t *wedge_new(int width, int height);

void wedge_free(wedge_t *w);

void wedge_processing_init(void);

void wedge_processing_queue(wedge_t *w);

int wedge_processing_dequeue(wedge_t *w);

#endif
