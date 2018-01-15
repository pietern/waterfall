#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <cairo/cairo.h>

#include "wedge.h"
#include "color.h"

// Only the address of this user data key is used.
// Contents is not relevant.
static cairo_user_data_key_t _data_key;

wedge_t *wedge_new(int width, int height) {
    wedge_t *w;

    w = calloc(1, sizeof(wedge_t));
    if (w == NULL) {
        assert(0);
    }

    w->width = width;
    w->height = height;

    return w;
}

static void wedge__willneed(wedge_t *w) {
    int64_t i;

    for (i = 0; i < w->height; i++) {
        stream_fadvise_willneed(w->stream, w->start + i * w->skip, w->width);
    }
}

static void wedge__dontneed(wedge_t *w) {
    int64_t i;

    for (i = 0; i < w->height; i++) {
        stream_fadvise_dontneed(w->stream, w->start + i * w->skip, w->width);
    }
}

void wedge_free(wedge_t *w) {
    wedge__dontneed(w);
    cairo_surface_destroy(w->surface);
    free(w);
}

static unsigned char *wedge__surface_alloc(wedge_t *w) {
    unsigned char *data;
    cairo_t *cr;

    w->stride = cairo_format_stride_for_width(CAIRO_FORMAT_ARGB32, w->width);

    data = calloc(w->stride, w->height);
    if (data == NULL) {
        assert(0);
    }

    w->surface = cairo_image_surface_create_for_data(data,
                                                     CAIRO_FORMAT_ARGB32,
                                                     w->width,
                                                     w->height,
                                                     w->stride);

    // Underlying buffer should be freed when surface is destroyed
    cairo_surface_set_user_data(w->surface, &_data_key, data, free);

    // Make surface transparent
    cr = cairo_create(w->surface);
    cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 0.0);
    cairo_paint(cr);
    cairo_destroy(cr);

    return data;
}

static void wedge_process(wedge_t *w) {
    stream_t *stream;
    double *mag;
    unsigned char *data;
    int64_t i;
    int64_t j;

    stream = w->stream;
    mag = malloc(sizeof(double) * w->width);
    data = wedge__surface_alloc(w);

    for (i = 0; i < w->height; i++) {
        if (w->start + i * w->skip + w->width >= stream->n_samples) {
            break;
        }

        memset(mag, 0, sizeof(double) * w->width);
        stream_fft(stream, w->start + i * w->skip, mag, w->width);

        // Normalize and populate wedge buffer
        int min = -100;
        int max = -60;
        for (j = 0; j < w->width; j++) {
            int16_t val;

            val = (255 * (mag[j] - min)) / (max - min);

            // Clamp
            if (val > 255) {
                val = 255;
            } else if (val < 0){
                val = 0;
            }

            data[i * w->stride + j*4 + 0] = color[val][2];
            data[i * w->stride + j*4 + 1] = color[val][1];
            data[i * w->stride + j*4 + 2] = color[val][0];
            data[i * w->stride + j*4 + 3] = 255;
        }
    }

    w->status = WEDGE_DONE;
    w->callback_fn(w, w->callback_data);

    free(mag);
}

static GList *queue = NULL;
static pthread_mutex_t lock;
static pthread_cond_t cond;

static void *wedge_worker(__attribute__((unused)) void *data) {
    GList *link;
    wedge_t *w;

    for (;;) {
        pthread_mutex_lock(&lock);
        for (;;) {
            link = g_list_first(queue);
            if (link == NULL) {
                pthread_cond_wait(&cond, &lock);
                continue;
            }

            w = (wedge_t *) link->data;
            queue = g_list_delete_link(queue, link);
            break;
        }

        pthread_mutex_unlock(&lock);

        w->status = WEDGE_PROCESSING;

        // Do actual work while not holding the lock
        wedge__willneed(w);
        wedge_process(w);
        wedge__dontneed(w);
    }

    return NULL;
}

void wedge_processing_init(void) {
    pthread_t t;

    queue = NULL;

    pthread_mutex_init(&lock, NULL);
    pthread_cond_init(&cond, NULL);
    pthread_create(&t, NULL, &wedge_worker, NULL);
}

void wedge_processing_queue(wedge_t *w) {
    w->status = WEDGE_QUEUED;

    pthread_mutex_lock(&lock);
    queue = g_list_append(queue, w);
    pthread_cond_signal(&cond);
    pthread_mutex_unlock(&lock);
}

int wedge_processing_dequeue(wedge_t *w) {
    GList *link;

    pthread_mutex_lock(&lock);
    link = g_list_find(queue, w);
    if (link == NULL) {
        pthread_mutex_unlock(&lock);
        return -1;
    }

    queue = g_list_remove_link(queue, link);
    pthread_mutex_unlock(&lock);
    return 0;
}
