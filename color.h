#ifndef _COLOR_H
#define _COLOR_H

#include <stdint.h>

typedef uint8_t color_t[3];

typedef color_t colormap_t[256];

// Pointer to first color in chosen colormap
extern color_t *color;

// Available colormaps
extern colormap_t colormap_rainbow;
extern colormap_t colormap_blue_to_yellow;

#endif
