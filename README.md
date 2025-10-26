# Forest Fire Simulation: Parallelism and Performance

This project contains three C programs for simulating a forest fire using a cellular automaton model, with core simulation logic parallelized using **OpenMP**.

## Project Files

* `forest_fire.c`: Graphical version (SDL3/OpenMP) with thread-local RNG.
* `forest_fire_time.c`: Performance timing version (OpenMP timers) with deterministic RNG for comparison.
* `forest_fire_serial.c` (Assumed): Serial baseline version.

***

## Compilation

### 1. Graphical Version (`forest_fire.c`)

Requires SDL3 and SDL3_image libraries.

```bash
gcc forest_fire.c -o forest_fire -I. -L. -lSDL3 -lSDL3_image -lm -fopenmp
```

### 2. Timing Version (`forest_fire_time.c`)

Requires OpenMP.

```bash
gcc forest_fire_time.c -o forest_fire_time -lm -fopenmp
```

-----

## Execution

The general command format is: `<executable> <grid_size> <wind_strength(0-1)> <wind_deg> [fire_x1,fire_y1 ...]`

### 1. Graphical Simulation (`forest_fire`)

  * **Linux/Unix Example:** `./forest_fire 100 0.9 0.0 50,50 20,20`
  * **Windows CMD Example:** `forest_fire.exe 100 0.9 0.0 50,50 20,20`

### 2. Performance Timing (`forest_fire_time`)

Can optionally take the number of threads as the last argument.

  * **Linux/Unix Example:** `./forest_fire_time 1000 0.5 90.0 4`
  * **Windows CMD Example:** `forest_fire_time.exe 1000 0.5 90.0 4`

-----

## Student Contributors

  * Hima Prasobh (2023BCS0086)
  * Rinu Ann Varghese (2023BCS0029)
