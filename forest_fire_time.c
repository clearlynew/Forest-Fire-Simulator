// forest_fire_time_deterministic.c
// SDL3 + SDL3_image + OpenMP forest fire simulation
// Deterministic version: fixed RNG based on cell position and step number

#define SDL_DISABLE_OLDNAMES 1
#include "SDL3/SDL.h"
#include "SDL3_image/SDL_image.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#else
static double omp_get_wtime(void){ return (double)clock() / CLOCKS_PER_SEC; }
static int omp_get_max_threads(void){ return 1; }
static int omp_get_thread_num(void){ return 0; }
static void omp_set_num_threads(int n) { (void)n; }
#endif

typedef enum { CELL_TREE = 0, CELL_BURNING = 1, CELL_BURNT = 2 } CellState;

// Deterministic RNG per cell and timestep
static inline double cell_rng(int x, int y, int step, unsigned int base_seed) {
    unsigned int seed = base_seed;
    seed ^= (unsigned int)(x * 73856093u);
    seed ^= (unsigned int)(y * 19349663u);
    seed ^= (unsigned int)(step * 83492791u);
    seed = (seed * 1103515245u + 12345u) & 0x7fffffffu;
    return (double)seed / 2147483647.0;
}

typedef struct {
    int grid_size;
    float wind_strength;
    double wind_deg;
    double wind_x, wind_y;
    CellState *cells;
    CellState *cells_next;
    unsigned int base_seed; // fixed deterministic seed
} Simulation;

int sim_init(Simulation *s, int grid_size, float wind_strength, double wind_deg,
             unsigned int base_seed, int *fire_coords, int fire_count) {
    s->grid_size = grid_size;
    s->wind_strength = fmaxf(0.0f, fminf(1.0f, wind_strength));
    s->wind_deg = wind_deg;
    double rad = wind_deg * M_PI / 180.0;
    s->wind_x = cos(rad);
    s->wind_y = sin(rad);
    s->base_seed = base_seed;

    size_t N = (size_t)grid_size * (size_t)grid_size;
    s->cells = (CellState*)malloc(N * sizeof(CellState));
    s->cells_next = (CellState*)malloc(N * sizeof(CellState));
    if (!s->cells || !s->cells_next) return 0;

    for (size_t i = 0; i < N; ++i)
        s->cells[i] = CELL_TREE;

    if (fire_count > 0 && fire_coords) {
        for (int p = 0; p < fire_count; ++p) {
            int x = fire_coords[p*2], y = fire_coords[p*2+1];
            if (x >= 0 && x < grid_size && y >= 0 && y < grid_size)
                s->cells[y * grid_size + x] = CELL_BURNING;
        }
    } else {
        s->cells[(grid_size/2) * grid_size + (grid_size/2)] = CELL_BURNING;
    }

    return 1;
}

void sim_free(Simulation *s) {
    if (s->cells) free(s->cells);
    if (s->cells_next) free(s->cells_next);
}

// One step of the deterministic fire spread
int sim_step(Simulation *s, int parallel, int step_num) {
    int G = s->grid_size;
    const int dr[8] = {-1, 1, 0, 0, -1, -1, 1, 1};
    const int dc[8] = { 0, 0,-1, 1, -1,  1,-1, 1};
    int burning_count = 0;

    if (parallel) {
        #pragma omp parallel for collapse(2) schedule(static) reduction(+:burning_count)
        for (int y = 0; y < G; ++y) {
            for (int x = 0; x < G; ++x) {
                int idx = y * G + x;
                CellState cur = s->cells[idx];
                CellState next = cur;

                if (cur == CELL_BURNING) {
                    next = CELL_BURNT;
                    burning_count++;
                } else if (cur == CELL_TREE) {
                    int ignite = 0;
                    for (int k = 0; k < 8; ++k) {
                        int ny = y + dr[k], nx = x + dc[k];
                        if ((unsigned)ny >= (unsigned)G || (unsigned)nx >= (unsigned)G)
                            continue;
                        if (s->cells[ny * G + nx] != CELL_BURNING)
                            continue;

                        double vx = x - nx, vy = y - ny;
                        double dist = sqrt(vx*vx + vy*vy);
                        if (dist <= 0.0) continue;
                        vx /= dist; vy /= dist;
                        double dot = vx * s->wind_x + vy * s->wind_y;
                        double base_p = 0.28;
                        double p = base_p * (1.0 + s->wind_strength * dot);
                        p = fmax(0.0, fmin(1.0, p));

                        // Deterministic per-cell RNG
                        if (cell_rng(x, y, step_num, s->base_seed) < p) {
                            ignite = 1;
                            break;
                        }
                    }
                    if (ignite) next = CELL_BURNING;
                }
                s->cells_next[idx] = next;
            }
        }
    } else { // serial
        for (int y = 0; y < G; ++y) {
            for (int x = 0; x < G; ++x) {
                int idx = y * G + x;
                CellState cur = s->cells[idx];
                CellState next = cur;

                if (cur == CELL_BURNING) {
                    next = CELL_BURNT;
                    burning_count++;
                } else if (cur == CELL_TREE) {
                    int ignite = 0;
                    for (int k = 0; k < 8; ++k) {
                        int ny = y + dr[k], nx = x + dc[k];
                        if ((unsigned)ny >= (unsigned)G || (unsigned)nx >= (unsigned)G)
                            continue;
                        if (s->cells[ny * G + nx] != CELL_BURNING)
                            continue;

                        double vx = x - nx, vy = y - ny;
                        double dist = sqrt(vx*vx + vy*vy);
                        if (dist <= 0.0) continue;
                        vx /= dist; vy /= dist;
                        double dot = vx * s->wind_x + vy * s->wind_y;
                        double base_p = 0.28;
                        double p = base_p * (1.0 + s->wind_strength * dot);
                        p = fmax(0.0, fmin(1.0, p));

                        //  Deterministic RNG again
                        if (cell_rng(x, y, step_num, s->base_seed) < p) {
                            ignite = 1;
                            break;
                        }
                    }
                    if (ignite) next = CELL_BURNING;
                }
                s->cells_next[idx] = next;
            }
        }
    }

    CellState *tmp = s->cells;
    s->cells = s->cells_next;
    s->cells_next = tmp;
    return burning_count;
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <grid_size> <wind_strength> <wind_deg> [num_threads]\n", argv[0]);
        return 1;
    }

    int grid_size = atoi(argv[1]);
    if (grid_size < 8) grid_size = 8;
    if (grid_size > 1000) grid_size = 1000;

    float wind_strength = (float)atof(argv[2]);
    double wind_deg = atof(argv[3]);

    int num_threads = 0;
    if (argc >= 5) {
        num_threads = atoi(argv[4]);
        if (num_threads > 0)
            omp_set_num_threads(num_threads);
    }

    printf("Using %d thread(s) for parallel simulation.\n", omp_get_max_threads());

    Simulation sim_serial, sim_parallel;
    int fire_point[2] = {grid_size/2, grid_size/2};

    unsigned int fixed_seed = 123456789u; //  constant seed

    sim_init(&sim_serial, grid_size, wind_strength, wind_deg, fixed_seed, fire_point, 1);
    sim_init(&sim_parallel, grid_size, wind_strength, wind_deg, fixed_seed, fire_point, 1);

    // --- Serial timing ---
    double t0 = omp_get_wtime();
    int step_s = 0, burning_s;
    do {
        burning_s = sim_step(&sim_serial, 0, step_s);
        step_s++;
    } while (burning_s > 0);
    double t1 = omp_get_wtime();
    printf("Serial finished in %d steps, time: %.6f s\n", step_s, t1 - t0);

    // --- Parallel timing ---
    double t2 = omp_get_wtime();
    int step_p = 0, burning_p;
    do {
        burning_p = sim_step(&sim_parallel, 1, step_p);
        step_p++;
    } while (burning_p > 0);
    double t3 = omp_get_wtime();
    printf("Parallel finished in %d steps, time: %.6f s\n", step_p, t3 - t2);

    printf("Hima Prasobh 2023BCS0086 \n Rinu Ann Varghese 2023BCS0029 ");

    sim_free(&sim_serial);
    sim_free(&sim_parallel);
    return 0;
}
