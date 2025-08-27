// forest_fire.c
// SDL3 + SDL3_image + OpenMP forest fire simulation
// Now with support for arbitrary/multiple fire starting points (x,y...) via command line

#define SDL_DISABLE_OLDNAMES 1  // Disable old SDL3 macros
#include "SDL3/SDL.h"
#include "SDL3_image/SDL_image.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#else
static int omp_get_max_threads(void){ return 1; }
static int omp_get_thread_num(void){ return 0; }
#endif

#define WIN_W 800
#define WIN_H 800

typedef enum { CELL_TREE = 0, CELL_BURNING = 1, CELL_BURNT = 2 } CellState;

/* Thread-local LCG RNG (31-bit) */
static inline unsigned int lcg_next(unsigned int *seed) {
    *seed = (*seed * 1103515245u + 12345u) & 0x7fffffffu;
    return *seed;
}
static inline double lcg_uniform(unsigned int *seed) {
    return (double)lcg_next(seed) / 2147483647.0;
}

/* Simulation state */
typedef struct {
    int grid_size;
    float wind_strength;
    double wind_deg;
    double wind_x, wind_y;
    CellState *cells;
    CellState *cells_next;
    unsigned int *thread_seeds;
} Simulation;

/* Load arrow PNG path */
static const char *ARROW_PATH = "blue-arrow.png";

/* Initialize simulation */
int sim_init(Simulation *s, int grid_size, float wind_strength, double wind_deg,
             unsigned int global_seed, int *fire_coords, int fire_count) {
    s->grid_size = grid_size;
    s->wind_strength = fmaxf(0.0f, fminf(1.0f, wind_strength));
    s->wind_deg = wind_deg;
    double rad = wind_deg * M_PI / 180.0;
    s->wind_x = cos(rad);
    s->wind_y = sin(rad);
    size_t N = (size_t)grid_size * (size_t)grid_size;
    s->cells = (CellState*)malloc(N * sizeof(CellState));
    s->cells_next = (CellState*)malloc(N * sizeof(CellState));
    if (!s->cells || !s->cells_next) return 0;
    for (size_t i = 0; i < N; ++i) s->cells[i] = CELL_TREE;
    // Ignite user-specified points, if any
    if (fire_count > 0 && fire_coords) {
        for (int p = 0; p < fire_count; ++p) {
            int x = fire_coords[p*2], y = fire_coords[p*2+1];
            if (x >= 0 && x < grid_size && y >= 0 && y < grid_size)
                s->cells[y * grid_size + x] = CELL_BURNING;
        }
    } else {
        // Default: ignite center
        s->cells[(grid_size/2) * grid_size + (grid_size/2)] = CELL_BURNING;
    }
    int T = omp_get_max_threads();
    s->thread_seeds = (unsigned int*)malloc(sizeof(unsigned int) * T);
    unsigned int base = global_seed ? global_seed : (unsigned int)time(NULL);
    for (int t = 0; t < T; ++t) s->thread_seeds[t] = base + 1337u * (t+1);
    return 1;
}
void sim_free(Simulation *s) {
    if (s->cells) free(s->cells);
    if (s->cells_next) free(s->cells_next);
    if (s->thread_seeds) free(s->thread_seeds);
}
/* Update one simulation step */
int sim_step(Simulation *s) {
    int G = s->grid_size;
    const int dr[8] = {-1, 1, 0, 0, -1, -1, 1, 1};
    const int dc[8] = { 0, 0,-1, 1, -1,  1,-1, 1};
    int burning_count = 0;
    #pragma omp parallel for collapse(2) schedule(static) reduction(+:burning_count)
    for (int y = 0; y < G; ++y) {
        for (int x = 0; x < G; ++x) {
            int idx = y * G + x;
            CellState cur = s->cells[idx];
            CellState next = cur;
            int tid = omp_get_thread_num();
            unsigned int seed = s->thread_seeds[tid];
            if (cur == CELL_BURNING) {
                next = CELL_BURNT;
                burning_count++;
            } else if (cur == CELL_TREE) {
                int ignite = 0;
                for (int k = 0; k < 8; ++k) {
                    int ny = y + dr[k], nx = x + dc[k];
                    if ((unsigned)ny >= (unsigned)G || (unsigned)nx >= (unsigned)G) continue;
                    if (s->cells[ny*G + nx] != CELL_BURNING) continue;
                    double vx = x - nx;
                    double vy = y - ny;
                    double dist = sqrt(vx*vx + vy*vy);
                    if (dist <= 0.0) continue;
                    vx /= dist; vy /= dist;
                    double dot = vx * s->wind_x + vy * s->wind_y;
                    double base_p = 0.28;
                    double p = base_p * (1.0 + s->wind_strength * dot);
                    p = fmax(0.0, fmin(1.0, p));
                    if (lcg_uniform(&seed) < p) { ignite = 1; break; }
                }
                if (ignite) next = CELL_BURNING;
            }
            s->cells_next[idx] = next;
            s->thread_seeds[tid] = seed;
        }
    }
    CellState *tmp = s->cells; s->cells = s->cells_next; s->cells_next = tmp;
    return burning_count;
}
/* Render everything */
void render_all(SDL_Renderer *renderer, Simulation *s, SDL_Texture *arrow_tex) {
    int G = s->grid_size;
    float cell_size = (float)WIN_W / G;
    SDL_SetRenderDrawColor(renderer, 240, 240, 240, 255);
    SDL_RenderClear(renderer);
    for (int y = 0; y < G; ++y) {
        for (int x = 0; x < G; ++x) {
            CellState st = s->cells[y*G + x];
            if (st == CELL_TREE) SDL_SetRenderDrawColor(renderer, 20, 150, 20, 255);
            else if (st == CELL_BURNING) SDL_SetRenderDrawColor(renderer, 255, 140, 0, 255);
            else SDL_SetRenderDrawColor(renderer, 30, 30, 30, 255);
            SDL_FRect r = { (float)x * cell_size, (float)y * cell_size, cell_size, cell_size };
            SDL_RenderFillRect(renderer, &r);
            SDL_SetRenderDrawColor(renderer, 0,0,0,80);
            SDL_RenderRect(renderer, &r);
        }
    }
    // Legend boxes
    int lx = WIN_W - 180;
    float ly = 20.0f;
    SDL_FRect box = { (float)lx, ly, 18.0f, 18.0f };
    SDL_SetRenderDrawColor(renderer, 20,150,20,255); SDL_RenderFillRect(renderer,&box);
    SDL_SetRenderDrawColor(renderer,0,0,0,255); SDL_RenderRect(renderer,&box);
    box.y += 26.0f;
    SDL_SetRenderDrawColor(renderer,255,140,0,255); SDL_RenderFillRect(renderer,&box);
    SDL_SetRenderDrawColor(renderer,0,0,0,255); SDL_RenderRect(renderer,&box);
    box.y += 26.0f;
    SDL_SetRenderDrawColor(renderer,30,30,30,255); SDL_RenderFillRect(renderer,&box);
    SDL_SetRenderDrawColor(renderer,255,255,255,255); SDL_RenderRect(renderer,&box);
    if (arrow_tex) {
        SDL_FRect dst = { (float)(WIN_W-120), 20.0f, 96.0f, 96.0f };
        double angle = s->wind_deg;
        SDL_RenderTextureRotated(renderer, arrow_tex, NULL, &dst, angle, NULL, SDL_FLIP_NONE);
    }
    SDL_RenderPresent(renderer);
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <grid_size> <wind_strength(0-1)> <wind_deg> [x1,y1 x2,y2 ...]\n", argv[0]);
        return 1;
    }

    // Corrected: use argv[1], argv[2], argv[3]
    int grid_size = atoi(argv[1]);
    if (grid_size < 8) grid_size = 8;
    if (grid_size > 1000) grid_size = 1000;

    float wind_strength = (float)atof(argv[2]);
    double wind_deg = atof(argv[3]);

    // Parse fire start points
    int fire_count = argc - 4;
    int *fire_coords = NULL;
    if (fire_count > 0) {
        fire_coords = malloc(fire_count * 2 * sizeof(int));
        for (int i = 0; i < fire_count; ++i) {
            int x = -1, y = -1;
            if (sscanf(argv[4 + i], "%d,%d", &x, &y) == 2
                && x >= 0 && x < grid_size && y >= 0 && y < grid_size) {
                fire_coords[i*2] = x;
                fire_coords[i*2+1] = y;
            } else {
                fprintf(stderr, "Invalid fire point: %s\n", argv[4 + i]);
                fire_coords[i*2] = fire_coords[i*2+1] = -1;
            }
        }
    }

    Simulation sim;
    if (!sim_init(&sim, grid_size, wind_strength, wind_deg, (unsigned int)time(NULL), fire_coords, fire_count)) {
        fprintf(stderr, "Failed to initialize simulation\n");
        if (fire_coords) free(fire_coords);
        return 1;
    }

    // SDL3: Returns true on success, false on failure!
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        fprintf(stderr, "SDL_Init failed: %s\n", SDL_GetError());
        sim_free(&sim);
        if (fire_coords) free(fire_coords);
        return 1;
    }
    SDL_Window *win = SDL_CreateWindow("Forest Fire (parallel)", WIN_W, WIN_H, SDL_WINDOW_RESIZABLE);
    if (!win) { fprintf(stderr,"SDL_CreateWindow failed: %s\n", SDL_GetError()); SDL_Quit(); sim_free(&sim); if (fire_coords) free(fire_coords); return 1; }
    SDL_Renderer *rend = SDL_CreateRenderer(win, NULL);
    if (!rend) { fprintf(stderr,"SDL_CreateRenderer failed: %s\n", SDL_GetError()); SDL_DestroyWindow(win); SDL_Quit(); sim_free(&sim); if (fire_coords) free(fire_coords); return 1; }
    SDL_Texture *arrow_tex = IMG_LoadTexture(rend, ARROW_PATH);
    if (!arrow_tex) {
        fprintf(stderr,"Warning: could not load arrow PNG at '%s'\n", ARROW_PATH);
    }

    int running = 1;
    SDL_Event e;
    int step = 0;
    const int ms_per_frame = 1000 / 15;
    while (running) {
        Uint64 t0 = SDL_GetTicks();
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_EVENT_QUIT || e.type == SDL_EVENT_WINDOW_CLOSE_REQUESTED) running = 0;
        }
        int burning = sim_step(&sim);
        render_all(rend, &sim, arrow_tex);
        step++;
        if (burning == 0) {
            printf("All burning cells extinguished at step %d.\n", step);
            SDL_Delay(800);
            break;
        }
        Uint64 t1 = SDL_GetTicks();
        int elapsed = (int)(t1 - t0);
        if (elapsed < ms_per_frame) SDL_Delay(ms_per_frame - elapsed);
    }
    if (arrow_tex) SDL_DestroyTexture(arrow_tex);
    SDL_DestroyRenderer(rend);
    SDL_DestroyWindow(win);
    SDL_Quit();
    sim_free(&sim);
    if (fire_coords) free(fire_coords);
    return 0;
}
