#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "window.h"
#include "stream.h"
#include "wedge.h"
#include "common.h"

#define SCALE_WIDTH 30
#define TICK_MARGIN 3

// User coordinates of viewport
typedef struct {
    int32_t w;
    int32_t h;

    // User coordinates of (0, 0) and (w, h)
    double uox;
    double uoy;
    double udx;
    double udy;
    double uw;
    double uh;
    double updx;
    double updy;

    // Scaled coordinates of (0, 0) and (w, h)
    double sox;
    double soy;
    double sdx;
    double sdy;
    double sw;
    double sh;
    double spdx;
    double spdy;

    // Scale factor in x and y axis
    int8_t sx;
    int8_t sy;

    // Offset of first wedge IN view
    int64_t start;

    // Offset of first wedge OUT of view
    int64_t end;

    // Separation between wedges
    int64_t skip;

    // Offsets of stream
    int64_t stream_start;
    int64_t stream_end;
} viewport_t;

typedef struct {
    uint16_t wedge_height;

    // Used to track dragging motion
    int xstart;
    int ystart;

    cairo_matrix_t window_to_user;
    cairo_matrix_t window_to_user_inv;

    cairo_matrix_t user_to_scaled;
    cairo_matrix_t user_to_scaled_inv;
    cairo_matrix_t user_to_scaled_start;

    viewport_t viewport;

    stream_t *stream;
    GPtrArray *wedges;

    uint64_t draw_cb_calls;

    GtkWidget *drawing_area;
    GtkWidget *time_scale;
    GtkWidget *frequency_scale;
} window_state_t;

static void _cairo_matrix_init_inverse(cairo_matrix_t *dst, const cairo_matrix_t *src) {
    memcpy(dst, src, sizeof(cairo_matrix_t));
    cairo_matrix_invert(dst);
}

static uint64_t timediff(struct timeval t1, struct timeval t2) {
    uint64_t m1 = t1.tv_sec * 1000000 + t1.tv_usec;
    uint64_t m2 = t2.tv_sec * 1000000 + t2.tv_usec;
    return m2 - m1;
}

static void compute_scale_factor(window_state_t *state, int8_t *sx, int8_t *sy) {
    double ddx = 1.0f;
    double ddy = 1.0f;
    double dfx;
    double dfy;

    // Find equivalent for 1x1 pixel block in user coordinates
    cairo_matrix_transform_distance(&state->window_to_user_inv, &ddx, &ddy);
    cairo_matrix_transform_distance(&state->user_to_scaled_inv, &ddx, &ddy);

    ddx = log2f(fabs(ddx));
    dfx = -floor(ddx);
    if (dfx < 0) {
        *sx = 0;
    } else {
        *sx = dfx;
    }

    ddy = log2f(fabs(ddy));
    dfy = -floor(ddy);
    if (dfy < 0) {
        *sy = 0;
    } else {
        *sy = dfy;
    }

    // At scale 14 the wedge width is 16384, which is the maximum cairo surface width.
    // (the cairo status flips to invalid value after using one of higher width)
    if (*sx > 14) {
        *sx = 14;
    }

    if (*sy > 14) {
        *sy = 14;
    }

    if (DEBUG) {
        fprintf(stderr, "sx: %d (ddx: %.3f), sy: %d (ddy: %.3f)\n", *sx, ddx, *sy, ddy);
    }

    return;
}

static void compute_viewport(window_state_t *state, viewport_t *v) {
    stream_t *stream = state->stream;
    double min;
    double max;
    double unit;

    compute_scale_factor(state, &v->sx, &v->sy);

    v->w = gtk_widget_get_allocated_width(state->drawing_area);
    v->h = gtk_widget_get_allocated_height(state->drawing_area);

    v->uox = 0.0f;
    v->uoy = 0.0f;
    v->udx = v->w;
    v->udy = v->h;
    v->uw = v->w;
    v->uh = v->h;
    v->updx = 1.0f;
    v->updy = 1.0f;

    // Find equivalent for the origin and window size in user coordinates
    cairo_matrix_transform_point(&state->window_to_user_inv, &v->uox, &v->uoy);
    cairo_matrix_transform_point(&state->window_to_user_inv, &v->udx, &v->udy);
    cairo_matrix_transform_distance(&state->window_to_user_inv, &v->uw, &v->uh);
    cairo_matrix_transform_distance(&state->window_to_user_inv, &v->updx, &v->updy);

    v->sox = v->uox;
    v->soy = v->uoy;
    v->sdx = v->udx;
    v->sdy = v->udy;
    v->sw = v->uw;
    v->sh = v->uh;
    v->spdx = v->updx;
    v->spdy = v->updy;

    // Find equivalent for the origin and window size in scaled coordinates
    cairo_matrix_transform_point(&state->user_to_scaled_inv, &v->sox, &v->soy);
    cairo_matrix_transform_point(&state->user_to_scaled_inv, &v->sdx, &v->sdy);
    cairo_matrix_transform_distance(&state->user_to_scaled_inv, &v->sw, &v->sh);
    cairo_matrix_transform_distance(&state->user_to_scaled_inv, &v->spdx, &v->spdy);

    if (v->soy < v->sdy) {
        min = v->soy;
        max = v->sdy;
    } else {
        min = v->sdy;
        max = v->soy;
    }

    // Determine offsets of first wedge in view and first wedge out of view
    unit = state->wedge_height * stream->sample_rate / powf(2.0f, v->sy);
    v->start = floor((min * stream->sample_rate) / unit) * unit;
    v->end = (floor((max * stream->sample_rate) / unit) + 1.0f) * unit;
    v->skip = unit;

    // Offsets of stream
    v->stream_start = 0;
    v->stream_end = stream->n_samples;
}

static void _wedge_row_callback(__attribute__((unused)) wedge_t *w, void *data) {
    GtkWidget *widget = (GtkWidget *) data;

    gtk_widget_queue_draw(widget);
}

static wedge_t *create_wedge(GtkWidget *widget, window_state_t *state, viewport_t *v, int64_t start, int64_t end) {
    wedge_t *w;

    w = wedge_new(1 << v->sx, state->wedge_height);

    w->stream = state->stream;

    w->sx = v->sx;
    w->sy = v->sy;
    w->start = start;
    w->end = end;
    w->skip = (end - start) / w->height;

    w->callback_fn = &_wedge_row_callback;
    w->callback_data = (void *) widget;

    wedge_processing_queue(w);

    if (DEBUG) {
        fprintf(stderr, "create wedge: sx: %d, sy: %d, start: %ld, end: %ld\n", w->sx, w->sy, w->start, w->end);
    }

    return w;
}

static void create_wedges_in_view(GtkWidget *widget, window_state_t *state) {
    viewport_t *v;
    stream_t *stream;
    int64_t i;
    int64_t j;
    wedge_t *w;

    compute_viewport(state, &state->viewport);

    v = &state->viewport;
    stream = state->stream;

    // Remove wedges that are queued for processing and are no longer
    // at the current scale level (e.g. when user rapidly zooms in/out).
    for (i = 0; i < state->wedges->len; i++) {
        w = (wedge_t *) g_ptr_array_index(state->wedges, i);
        if (w->status != WEDGE_QUEUED) {
            continue;
        }

        // Skip queued wedges at the right scale level
        if (w->sx == v->sx && w->sy == v->sy) {
            continue;
        }

        // Attempt to dequeue wedge for processing.
        // This can race with the processing thread so might fail.
        if (wedge_processing_dequeue(w) == 0) {
            g_ptr_array_remove_index(state->wedges, i--);
            wedge_free(w);
        }
    }

    // Remove wedges that were not used in the most recent draw call.
    for (i = 0; i < state->wedges->len; i++) {
        w = (wedge_t *) g_ptr_array_index(state->wedges, i);
        if (w->status != WEDGE_DONE) {
            continue;
        }

        // At the current scale level, never remove wedges that
        // are in view or are neighboring the view.
        if (w->sx == v->sx &&
            w->sy == v->sy &&
            w->start >= (v->start - v->skip) &&
            w->start < (v->end + v->skip))
        {
            continue;
        }

        // If it wasn't used in the most recent draw call, remove.
        if (w->draw_cb < state->draw_cb_calls) {
            g_ptr_array_remove_index(state->wedges, i--);
            wedge_free(w);
        }
    }

    // Iterate over wedges for this view.
    // Note that both start and end are extended by skip so that
    // wedges just outside the view port will be created before they
    // are visible, reducing the amount of flickering on screen.
    for (i = (v->start - v->skip); i < (v->end + v->skip); i += v->skip) {
        // Skip negative start offset
        if (i < 0) {
            continue;
        }

        // Break on out of bounds start offset
        if (i >= stream->n_samples) {
            break;
        }

        // Try to find the same wedge
        for (j = 0; j < state->wedges->len; j++) {
            w = g_ptr_array_index(state->wedges, j);
            if (w->sx == v->sx && w->sy == v->sy && w->start == i) {
                break;
            }
        }

        // Skip if this wedge was already computed
        if (j < state->wedges->len) {
            continue;
        }

        w = create_wedge(widget, state, v, i, i + v->skip);
        g_ptr_array_add(state->wedges, w);
    }

    return;
}

static void configure_event_cb (__attribute__((unused)) GtkWidget *widget,
                                __attribute__((unused)) GdkEventConfigure *event,
                                __attribute__((unused)) gpointer data) {
    window_state_t *state;
    static double sx = NAN;
    static double sy = NAN;
    double dx;
    double dy;
    int w = gtk_widget_get_allocated_width(widget);

    state = (window_state_t *) data;

    // Initialize sx and sy on the first configuration event.
    // Don't do this on subsequent events to avoid scaling on resize.
    if (isnan(sx) && isnan(sy)) {
        sx = w;
        sy = -1.0f;
    }

    dx = w / sx / 2.0f;
    dy = -state->stream->n_samples / state->stream->sample_rate;

    // Recompute window_to_user matrix
    cairo_matrix_init_identity(&state->window_to_user);
    cairo_matrix_scale(&state->window_to_user, sx, sy);
    cairo_matrix_translate(&state->window_to_user, dx, dy);

    // Recompute inverse
    _cairo_matrix_init_inverse(&state->window_to_user_inv, &state->window_to_user);

    // Ensure all wedges in view are available
    create_wedges_in_view(widget, state);

    gtk_widget_queue_draw(widget);
}

static gint _wedge_sort_fn(gconstpointer a, gconstpointer b) {
    wedge_t *wa = *((wedge_t **) a);
    wedge_t *wb = *((wedge_t **) b);

    // Check explicitly instead of returning 'wb - wa' to avoid
    // the gint from overflowing.
    if (wa->start < wb->start) {
        return -1;
    } else if (wa->start > wb->start) {
        return 1;
    } else {
        return 0;
    }
}

static gint _wedges_cover_viewport(GPtrArray *wedges, const viewport_t *v) {
    int64_t start = -1;
    int64_t end = -1;
    uint8_t i;
    wedge_t *w;
    int64_t start_max;
    int64_t end_min;

    g_ptr_array_sort(wedges, _wedge_sort_fn);

    // Check if wedges are consecutive
    for (i = 0; i < wedges->len; i++) {
        w = (wedge_t *) g_ptr_array_index(wedges, i);
        if (i == 0) {
            start = w->start;
            end = w->end;
        } else {
            if (w->start != end) {
                return -1;
            }

            end = w->end;
        }
    }

    // Check if range covers viewport
    start_max = v->stream_start;
    if (v->start > start_max) {
        start_max = v->start;
    }
    end_min = v->stream_end;
    if (v->end < end_min) {
        end_min = v->end;
    }
    if (start > start_max || end < end_min) {
        return -1;
    }

    return 0;
}

// Wedges are accumulated from best fit to worst fit (with respect to their scale).
// To draw them as a stack starting with the worst fit and ending with the
// best fit, the returned array must be traversed from the last element to the first.
static GPtrArray *_window_find_wedges_to_draw(window_state_t *state, const viewport_t *v) {
    GPtrArray *wedges_to_draw;
    GPtrArray *wedges_for_scale;
    int8_t sx;
    int8_t sy;
    int8_t dx;
    int8_t dy;
    int8_t wi;
    int8_t wn;
    uint16_t i;
    wedge_t *w;
    int rv;

    wedges_to_draw = g_ptr_array_new();

    // Draw wedges for scale level. If wedges at this level don't
    // cover the viewport, try progressively deviating scales.
    sx = v->sx;
    sy = v->sy;
    dx = 0;
    dy = 0;
    wi = 1;
    wn = 0;
    for (;;) {
        // Keep array of wedges at this scale level so we can
        // check later if the viewport was covered.
        wedges_for_scale = g_ptr_array_new();

        // Find applicable wedges
        for (i = 0; i < state->wedges->len; i++) {
            w = g_ptr_array_index(state->wedges, i);

            // Only draw processed wedges
            if (w->status != WEDGE_DONE) {
                continue;
            }

            // Scale must match
            if (w->sx != sx || w->sy != sy) {
                continue;
            }

            // Wedge must be in view
            if (w->end < v->start || w->start > v->end) {
                continue;
            }

            g_ptr_array_add(wedges_for_scale, w);
            g_ptr_array_add(wedges_to_draw, w);
        }

        // Check if drawn wedges cover the viewport
        rv = _wedges_cover_viewport(wedges_for_scale, v);
        g_ptr_array_unref(wedges_for_scale);
        if (rv == 0) {
            break;
        }

        // If distance from center scale (v.sx, v.sy) is done, grow range
        if (--wi <= 0) {
            // Stop after walking at distance 3
            if (wn++ >= 3) {
                break;
            }

            sy++;
            dx = 1;
            dy = 0;
            wi = 8 * wn;
            continue;
        }

        // Walk around (v->sx, v->sy) counter clockwise
        for (;;) {
            sx += dx;
            sy += dy;
            if (sy < v->sy - wn) {
                dx = 1;
                sx = sx + dx;
                dy = 0;
                sy = v->sy - wn;
            } else if (sy > v->sy + wn) {
                dx = -1;
                sx = sx + dx;
                dy = 0;
                sy = v->sy + wn;
            } else if (sx < v->sx - wn) {
                dx = 0;
                sx = v->sx - wn;
                dy = -1;
                sy = sy + dy;
            } else if (sx > v->sx + wn) {
                dx = 0;
                sx = v->sx + wn;
                dy = 1;
                sy = sy + dy;
            }

            if (sx >= 0 && sy >= 0) {
                break;
            } else {
                // Keep going if out of bounds
                continue;
            }
        }
    }

    return wedges_to_draw;
}

static gboolean draw_cb(__attribute__((unused)) GtkWidget *widget,
                        cairo_t *cr,
                        gpointer data) {
    window_state_t *state;
    GPtrArray *wedges;
    uint16_t i;
    wedge_t *w;

    state = (window_state_t *) data;
    wedges = _window_find_wedges_to_draw(state, &state->viewport);

    // Keep counter for draw callbacks
    state->draw_cb_calls++;

    // Apply transformations in accumulators and draw wedges
    cairo_transform(cr, &state->window_to_user);
    cairo_transform(cr, &state->user_to_scaled);

    // Paint background
    cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
    cairo_set_operator(cr, CAIRO_OPERATOR_OVER);
    cairo_paint(cr);

    // Paint patterns from base to current scale
    for (i = wedges->len; i > 0; i--) {
        double dx;
        double dy;

        w = g_ptr_array_index(wedges, i - 1);
        w->draw_cb = state->draw_cb_calls;

        dx = 0.0f;
        dy = (double) w->start / (double) state->stream->sample_rate;

        cairo_save(cr);

        cairo_translate(cr, dx, dy);
        cairo_scale(cr, 1.0f / w->width, 1.0f / powf(2.0f, w->sy));
        cairo_set_operator(cr, CAIRO_OPERATOR_OVER);
        cairo_set_source_surface(cr, w->surface, -(w->width / 2), 0);

        if (1) {
            cairo_pattern_set_extend(cairo_get_source(cr), CAIRO_EXTEND_PAD);
            cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_BILINEAR);
        } else {
            cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_NEAREST);
        }

        cairo_rectangle(cr, -(w->width / 2), -1, w->width, w->height + 1);
        cairo_fill(cr);

        cairo_restore(cr);
    }

    g_ptr_array_unref(wedges);

    return FALSE;
}

static void _to_human_time(char *buf, size_t len, double t, double dt) {
    int32_t h;
    int32_t m;
    int32_t s;
    int32_t ms;
    char min[2];
    uint8_t p = 0;

    // Round to nearest millisecond
    t = round(1000.0f * t) / 1000.0f;
    if (t < 0) {
        t *= -1.0f;
        min[0] = '-';
        min[1] = '\0';
    } else {
        min[0] = '\0';
    }

    // Compute fractions
    h = floor(t / 3600.0f);
    t -= h * 3600.0f;
    m = floor(t / 60.0f);
    t -= m * 60.0f;
    s = floor(t);
    t -= s;
    ms = round(t * 1000.0f);

    p += snprintf(buf + p, len - p, "%s%02d:%02d:%02d", min, h, m, s);
    if (dt < 1.0f) {
        p += snprintf(buf + p, len - p, ".%03d", ms);
    }
}

static void _window_time_scale_draw_axis(window_state_t *state, cairo_t *cr) {
    const viewport_t *v = &state->viewport;
    double scale = SCALE_WIDTH / 36.0f;
    double font_size = 11.0f * scale;

    // Use a minimum of 120 pixels between major ticks
    double dy = scale * 120.0f * fabs(v->spdy);
    double major_tick_delta;

    // Different staggering if > 1s or < 1s per major tick
    if (dy >= 1.0f) {
        double mul;

        // Make routine work for minutes and seconds
        mul = 1.0f;
        if (dy > 60.0f) {
            dy /= 60.0f;
            mul = 60.0f;
        }

        if (dy > 30.0f) {
            major_tick_delta = 60.0f * mul;
        } else if (dy > 10.0f) {
            major_tick_delta = 30.0f * mul;
        } else if (dy > 5.0f) {
            major_tick_delta = 10.0f * mul;
        } else if (dy > 2.0f) {
            major_tick_delta = 5.0f * mul;
        } else if (dy > 1.0f) {
            major_tick_delta = 2.0f * mul;
        } else {
            major_tick_delta = 1.0f * mul;
        }
    } else {
        double base;
        double rem;

        base = pow(10.0f, floor(log10(dy)));
        rem = dy / base;

        if (rem > 5.0f) {
            major_tick_delta = 10.0f * base;
        } else if (rem > 2.0f) {
            major_tick_delta = 5.0f * base;
        } else if (rem > 1.0f) {
            major_tick_delta = 2.0f * base;
        } else {
            major_tick_delta = 1.0f * base;
        }
    }

    // 10 minor ticks per major tick
    double minor_tick_delta = major_tick_delta / 10.0f;

    // The viewport is flipped on the x axis (increases from bottom to top)
    double t_min = v->sdy;
    double t_max = v->soy;

    // Start just outside the visible area so that text that overflows into
    // the visible area is still drawn.
    double tick = major_tick_delta * ceil(t_min / major_tick_delta) - major_tick_delta;

    // Continue into just outside the visible area.
    while (tick < (t_max + major_tick_delta)) {
        double x = 0.0f;
        double y = tick;
        char buf[32];

        _to_human_time(buf, sizeof(buf), tick, major_tick_delta);

        cairo_matrix_transform_point(&state->user_to_scaled, &x, &y);
        cairo_matrix_transform_point(&state->window_to_user, &x, &y);

        cairo_new_path(cr);
        cairo_move_to(cr, SCALE_WIDTH - TICK_MARGIN, y);
        cairo_line_to(cr, TICK_MARGIN, y);
        cairo_set_line_width(cr, scale);
        cairo_stroke(cr);

        cairo_save(cr);
        cairo_translate(cr, SCALE_WIDTH / 2.0f, y - 3);
        cairo_rotate(cr, -0.25 * M_TAU);
        cairo_set_font_size(cr, font_size);
        cairo_show_text(cr, buf);
        cairo_restore(cr);

        tick += major_tick_delta;
    }

    // Start just outside the visible area so that text that overflows into
    // the visible area is still drawn.
    tick = minor_tick_delta * (ceil(t_min / minor_tick_delta) - 1);

    // Continue into just outside the visible area.
    while (tick < (t_max + minor_tick_delta)) {
        double x = 0.0f;
        double y = tick;

        cairo_matrix_transform_point(&state->user_to_scaled, &x, &y);
        cairo_matrix_transform_point(&state->window_to_user, &x, &y);

        cairo_new_path(cr);
        cairo_move_to(cr, SCALE_WIDTH - TICK_MARGIN, y);
        cairo_line_to(cr, (2.0f * SCALE_WIDTH / 3.0f) - TICK_MARGIN, y);
        cairo_set_line_width(cr, scale);
        cairo_stroke(cr);

        tick += minor_tick_delta;
    }
}

static gboolean _window_time_scale_draw(__attribute__((unused)) GtkWidget *widget,
                                        cairo_t *cr,
                                        gpointer data) {
    window_state_t *state;
    GtkStyleContext *style;
    GdkRGBA bg;
    GdkRGBA fg;

    state = (window_state_t *) data;
    style = gtk_widget_get_style_context(widget);
    gtk_style_context_get_background_color(style, GTK_STATE_FLAG_NORMAL, &bg);
    gtk_style_context_get_color(style, GTK_STATE_FLAG_NORMAL, &fg);

    cairo_set_source_rgba(cr, bg.red, bg.green, bg.blue, bg.alpha);
    cairo_set_operator(cr, CAIRO_OPERATOR_OVER);
    cairo_paint(cr);

    cairo_set_source_rgba(cr, fg.red, fg.green, fg.blue, fg.alpha);
    _window_time_scale_draw_axis(state, cr);

    return FALSE;
}

static void _to_human_frequency(window_state_t *state, char *buf, size_t len, double t, double dt) {
    const char *suffix;
    double div;
    int8_t lcf = floor(log10(state->stream->center_freq));
    if (lcf >= 9) {
        suffix = "GHz";
        div = 1000000000.0f;
    } else if (lcf >= 6){
        suffix = "MHz";
        div = 1000000.0f;
    } else if (lcf >= 3) {
        suffix = "KHz";
        div = 1000.0f;
    } else {
        suffix = "Hz";
        div = 1.0f;
    }

    char fmt[16] = "%.0f %s";
    int8_t precision = (int8_t) (floor(log10(div)) - floor(log10(dt)));
    if (precision > 0) {
        fmt[2] = '0' + precision;
    }

    snprintf(buf, len, fmt, t / div, suffix);
}

static void _window_frequency_scale_draw_axis(window_state_t *state, cairo_t *cr) {
    const viewport_t *v = &state->viewport;
    double scale = SCALE_WIDTH / 36.0f;
    double font_size = 11.0f * scale;

    // Use a minimum of 120 pixels between major ticks
    double dx = scale * 120.0f * fabs(v->spdx);

    // Multiply by sample rate to normalize w.r.t. actual frequency,
    // not the fake scale of 1 unit per the complete spectrum.
    // The resulting scale is Hz between major ticks.
    dx *= state->stream->sample_rate;

    double base = pow(10.0f, floor(log10(dx)));
    double rem = dx / base;
    double major_tick_delta;

    if (rem > 5.0f) {
        major_tick_delta = 10.0f * base;
    } else if (rem > 2.0f) {
        major_tick_delta = 5.0f * base;
    } else if (rem > 1.0f) {
        major_tick_delta = 2.0f * base;
    } else {
        major_tick_delta = 1.0f * base;
    }

    // 10 minor ticks per major tick
    double minor_tick_delta = major_tick_delta / 10.0f;

    // Translate min/max to frequency in Hz.
    double t_min = state->stream->center_freq + state->stream->sample_rate * v->sox;
    double t_max = state->stream->center_freq + state->stream->sample_rate * v->sdx;

    // Start just outside the visible area so that text that overflows into
    double tick = major_tick_delta * ceil(t_min / major_tick_delta) - major_tick_delta;

    // Continue into just outside the visible area.
    while (tick < (t_max + major_tick_delta)) {
        double x = (tick - state->stream->center_freq) / state->stream->sample_rate;
        double y = 0.0f;
        char buf[32];

        _to_human_frequency(state, buf, sizeof(buf), tick, major_tick_delta);

        cairo_matrix_transform_point(&state->user_to_scaled, &x, &y);
        cairo_matrix_transform_point(&state->window_to_user, &x, &y);

        cairo_new_path(cr);
        cairo_move_to(cr, x, SCALE_WIDTH - TICK_MARGIN);
        cairo_line_to(cr, x, TICK_MARGIN);
        cairo_set_line_width(cr, scale);
        cairo_stroke(cr);

        cairo_save(cr);
        cairo_translate(cr, x + 3, SCALE_WIDTH / 2.0f);
        cairo_set_font_size(cr, font_size);
        cairo_show_text(cr, buf);
        cairo_restore(cr);

        tick += major_tick_delta;
    }

    // Start just outside the visible area so that text that overflows into
    // the visible area is still drawn.
    tick = minor_tick_delta * (ceil(t_min / minor_tick_delta) - 1);

    // Continue into just outside the visible area.
    while (tick < (t_max + minor_tick_delta)) {
        double x = (tick - state->stream->center_freq) / state->stream->sample_rate;
        double y = 0.0f;

        cairo_matrix_transform_point(&state->user_to_scaled, &x, &y);
        cairo_matrix_transform_point(&state->window_to_user, &x, &y);

        cairo_new_path(cr);
        cairo_move_to(cr, x, SCALE_WIDTH - TICK_MARGIN);
        cairo_line_to(cr, x, (2.0f * SCALE_WIDTH / 3.0f) - TICK_MARGIN);
        cairo_set_line_width(cr, scale);
        cairo_stroke(cr);

        tick += minor_tick_delta;
    }
}

static gboolean _window_frequency_scale_draw(__attribute__((unused)) GtkWidget *widget,
                                             cairo_t *cr,
                                             gpointer data) {
    window_state_t *state;
    GtkStyleContext *style;
    GdkRGBA bg;
    GdkRGBA fg;

    state = (window_state_t *) data;
    style = gtk_widget_get_style_context(widget);
    gtk_style_context_get_background_color(style, GTK_STATE_FLAG_NORMAL, &bg);
    gtk_style_context_get_color(style, GTK_STATE_FLAG_NORMAL, &fg);

    cairo_set_source_rgba(cr, bg.red, bg.green, bg.blue, bg.alpha);
    cairo_set_operator(cr, CAIRO_OPERATOR_OVER);
    cairo_paint(cr);

    cairo_set_source_rgba(cr, fg.red, fg.green, fg.blue, fg.alpha);
    _window_frequency_scale_draw_axis(state, cr);

    return FALSE;
}

static void close_window(void) {
    gtk_main_quit();
}

static void window_state_restore_start(window_state_t *state) {
    memcpy(&state->user_to_scaled, &state->user_to_scaled_start, sizeof(cairo_matrix_t));
    _cairo_matrix_init_inverse(&state->user_to_scaled_inv, &state->user_to_scaled);
}

static void window_state_translate(window_state_t *state,
                                   double dx,
                                   double dy)
{
    // Convert position delta to user coordinates to scaled coordinates
    cairo_matrix_transform_distance(&state->window_to_user_inv, &dx, &dy);
    cairo_matrix_transform_distance(&state->user_to_scaled_inv, &dx, &dy);

    // Apply transformation
    cairo_matrix_translate(&state->user_to_scaled, dx, dy);
    _cairo_matrix_init_inverse(&state->user_to_scaled_inv, &state->user_to_scaled);
}

static void window_state_scale(window_state_t *state,
                               double offset_x,
                               double offset_y,
                               double scale_x,
                               double scale_y)
{
    // Convert offset to user coordinates to scaled coordinates
    cairo_matrix_transform_point(&state->window_to_user_inv, &offset_x, &offset_y);
    cairo_matrix_transform_point(&state->user_to_scaled_inv, &offset_x, &offset_y);

    // Scale around pivot point equal to mouse position when mouse drag started
    cairo_matrix_translate(&state->user_to_scaled, offset_x, offset_y);
    cairo_matrix_scale(&state->user_to_scaled, scale_x, scale_y);
    cairo_matrix_translate(&state->user_to_scaled, -offset_x, -offset_y);
    _cairo_matrix_init_inverse(&state->user_to_scaled_inv, &state->user_to_scaled);
}

static void window_state_create_wedges(window_state_t *state,
                                       GtkWidget *widget)
{
    struct timeval t1;
    struct timeval t2;

    gettimeofday(&t1, NULL);
    create_wedges_in_view(widget, state);
    gettimeofday(&t2, NULL);

    gtk_widget_queue_draw(widget);
    gtk_widget_queue_draw(state->time_scale);
    gtk_widget_queue_draw(state->frequency_scale);

    if (DEBUG) {
        fprintf(stderr, "create_wedges: %luus\n", timediff(t1, t2));
    }
}

static gboolean motion_notify_event_cb (GtkWidget *widget,
                                        GdkEventMotion *event,
                                        gpointer data)
{
    window_state_t *state;

    state = (window_state_t *) data;

    if (event->state & GDK_BUTTON1_MASK) {
        double ddx = event->x - state->xstart;
        double ddy = event->y - state->ystart;

        window_state_restore_start(state);

        window_state_translate(state, ddx, ddy);

        window_state_create_wedges(state, widget);
    }

    if (event->state & GDK_BUTTON3_MASK) {
        int dx = event->x - state->xstart;
        int dy = event->y - state->ystart;
        double sdx = pow(2.0f, dx / 100.0f);
        double sdy = pow(2.0f, dy / 100.0f);
        double offset_x = state->xstart;
        double offset_y = state->ystart;

        window_state_restore_start(state);

        window_state_scale(state, offset_x, offset_y, sdx, sdy);

        window_state_create_wedges(state, widget);
    }

    return TRUE;
}

static gboolean button_press_event_cb (__attribute__((unused)) GtkWidget *widget,
                                       GdkEventButton *event,
                                       gpointer data)
{
    window_state_t *state;

    state = (window_state_t *) data;
    if (event->type == GDK_BUTTON_PRESS) {
        state->xstart = event->x;
        state->ystart = event->y;
        memcpy(&state->user_to_scaled_start, &state->user_to_scaled, sizeof(cairo_matrix_t));
    }

    return TRUE;
}

static gboolean key_press_event_cb(GtkWidget *widget,
                                   GdkEventKey *event,
                                   gpointer data)
{
    window_state_t *state;
    int shift;
    int w;
    int h;
    float scale = 2.0f;
    float translate = 25.0f;

    state = (window_state_t *) data;
    shift = event->state & GDK_SHIFT_MASK;
    w = gtk_widget_get_allocated_width(widget);
    h = gtk_widget_get_allocated_height(widget);

    if (shift) {
        // Shift modifier
        switch (event->keyval) {
        case GDK_KEY_Left:
            window_state_scale(state, w / 2, h / 2, 1.0f / scale, 1.0f);
            break;
        case GDK_KEY_Up:
            window_state_scale(state, w / 2, h / 2, 1.0f, 1.0f / scale);
            break;
        case GDK_KEY_Right:
            window_state_scale(state, w / 2, h / 2, scale, 1.0f);
            break;
        case GDK_KEY_Down:
            window_state_scale(state, w / 2, h / 2, 1.0f, scale);
            break;
        default:
            return FALSE;
        }

        window_state_create_wedges(state, widget);
    } else {
        // No modifier
        switch (event->keyval) {
        case GDK_KEY_Left:
            window_state_translate(state, translate, 0.0f);
            break;
        case GDK_KEY_Up:
            window_state_translate(state, 0.0f, translate);
            break;
        case GDK_KEY_Right:
            window_state_translate(state, -translate, 0.0f);
            break;
        case GDK_KEY_Down:
            window_state_translate(state, 0.0f, -translate);
            break;
        default:
            return FALSE;
        }

        window_state_create_wedges(state, widget);
    }

    return TRUE;
}

static GtkWidget *_create_drawing_areas(window_state_t *state) {
    GtkWidget *drawing_area;
    GtkWidget *time_scale;
    GtkWidget *frequency_scale;
    GtkWidget *grid;

    drawing_area = gtk_drawing_area_new();
    gtk_widget_set_hexpand(drawing_area, TRUE);
    gtk_widget_set_vexpand(drawing_area, TRUE);
    gtk_widget_set_size_request(drawing_area, 200, 200);

    // Signals used to handle the backing surface
    g_signal_connect(drawing_area, "draw",
                     G_CALLBACK(draw_cb), state);
    g_signal_connect(drawing_area, "configure-event",
                     G_CALLBACK(configure_event_cb), state);
    g_signal_connect(drawing_area, "motion-notify-event",
                     G_CALLBACK(motion_notify_event_cb), state);
    g_signal_connect(drawing_area, "button-press-event",
                     G_CALLBACK(button_press_event_cb), state);
    g_signal_connect(drawing_area, "key-press-event",
                     G_CALLBACK(key_press_event_cb), state);

    gtk_widget_set_events(drawing_area, gtk_widget_get_events(drawing_area)
                          | GDK_POINTER_MOTION_MASK
                          | GDK_BUTTON_PRESS_MASK
                          | GDK_KEY_PRESS_MASK);

    gtk_widget_set_can_focus(drawing_area, TRUE);

    time_scale = gtk_drawing_area_new();
    gtk_widget_set_hexpand(time_scale, FALSE);
    gtk_widget_set_vexpand(time_scale, TRUE);
    gtk_widget_set_size_request(time_scale, SCALE_WIDTH, SCALE_WIDTH);

    g_signal_connect(time_scale, "draw",
                     G_CALLBACK(_window_time_scale_draw), state);

    frequency_scale = gtk_drawing_area_new();
    gtk_widget_set_hexpand(frequency_scale, TRUE);
    gtk_widget_set_vexpand(frequency_scale, FALSE);
    gtk_widget_set_size_request(frequency_scale, SCALE_WIDTH, SCALE_WIDTH);

    g_signal_connect(frequency_scale, "draw",
                     G_CALLBACK(_window_frequency_scale_draw), state);

    grid = gtk_grid_new();
    gtk_grid_set_row_homogeneous(GTK_GRID(grid), FALSE);
    gtk_grid_set_column_homogeneous(GTK_GRID(grid), FALSE);

    gtk_grid_attach(GTK_GRID(grid), time_scale, 0, 1, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), frequency_scale, 1, 0, 1, 1);
    gtk_grid_attach(GTK_GRID(grid), drawing_area, 1, 1, 1, 1);

    state->drawing_area = drawing_area;
    state->time_scale = time_scale;
    state->frequency_scale = frequency_scale;

    return grid;
}

static void activate(GtkApplication *app,
                     gpointer user_data) {
    GtkWidget *window;
    GtkWidget *main_box;
    GtkWidget *drawing_area;
    stream_t *s;
    window_state_t *state;

    s = (stream_t *) user_data;

    state = malloc(sizeof(*state));
    state->wedge_height = 64;
    state->xstart = 0;
    state->ystart = 0;
    state->stream = s;
    state->wedges = g_ptr_array_new();

    cairo_matrix_init_identity(&state->user_to_scaled);
    _cairo_matrix_init_inverse(&state->user_to_scaled_inv, &state->user_to_scaled);

    window = gtk_application_window_new(app);
    gtk_window_set_title(GTK_WINDOW(window), "Waterfall");
    g_signal_connect(window, "destroy", G_CALLBACK(close_window), NULL);

    main_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
    gtk_container_add(GTK_CONTAINER(window), main_box);

    drawing_area = _create_drawing_areas(state);
    gtk_box_pack_start(GTK_BOX(main_box), drawing_area, TRUE, TRUE, 0);

    gtk_widget_show_all(window);
}

int window_run(stream_t *s) {
    GtkApplication *app;
    int status;

    app = gtk_application_new("com.github.pietern.waterfall", G_APPLICATION_FLAGS_NONE);
    g_signal_connect(app, "activate", G_CALLBACK(activate), s);
    status = g_application_run(G_APPLICATION(app), 0, NULL);
    g_object_unref(app);

    return status;
}
