#ifndef PERF_H_
#define PERF_H_

#include <swperf.h>
#include <string.h>

// extern void penv_slave0_cycle_init();
// extern void penv_slave0_cycle_count(unsigned long *count);

#define MAX_PERF_UNITS 32
#define MAX_CALL_STACK 32
#define PERF_INCLUSIVE 0
#define PERF_EXCLUSIVE 1

typedef struct {
    int slave_id;

    unsigned long clocks;
    int unit_num;
    const char *unit_name[MAX_PERF_UNITS];
    int type[MAX_PERF_UNITS];
    unsigned long total_clocks[MAX_PERF_UNITS];

    int call_stack[MAX_CALL_STACK];
    int stack_top;
} PerfEnv;

static inline void slave_init_perf_env(PerfEnv *env, int slave_id) {
    memset(env, 0, sizeof(PerfEnv));
    env->slave_id = slave_id;
    penv_slave0_cycle_init();
}

static inline int slave_perf_get_unit_id_by_name(PerfEnv *env, const char *name) {
    int unit_id = -1;
    for (int i = 0; i < env->unit_num; ++i) {
        if (strcmp(name, env->unit_name[i]) == 0) {
            unit_id = i;
            break;
        }
    }
    return unit_id;
}

static inline void slave_perf_checkpoint(PerfEnv *env) {
    unsigned long count;
    penv_slave0_cycle_count(&count);
    env->clocks += count;
    for (int i = 0; i < env->stack_top - 1; ++i) {
        int unit_id = env->call_stack[i];
        if (env->type[unit_id] == PERF_INCLUSIVE) {
            env->total_clocks[unit_id] += count;
        }
    }

    if (env->stack_top > 0) {
        int unit_id = env->call_stack[env->stack_top - 1];
        env->total_clocks[unit_id] += count;
    }

    penv_slave0_cycle_init();
}

static inline void slave_perf_begin(PerfEnv *env, int type, const char *name) {
    slave_perf_checkpoint(env);
    int unit_id = slave_perf_get_unit_id_by_name(env, name);
    if (unit_id == -1) {
        unit_id = env->unit_num++;
        env->unit_name[unit_id] = name;
        env->type[unit_id] = type;
    }
    env->call_stack[env->stack_top++] = unit_id;

    penv_slave0_cycle_init();
}

static inline void slave_perf_end(PerfEnv *env, const char *name) {
    slave_perf_checkpoint(env);
    int unit_id = slave_perf_get_unit_id_by_name(env, name);
    env->stack_top -= 1;

    penv_slave0_cycle_init();
}

static inline void slave_perf_report(PerfEnv *env) {
    slave_perf_checkpoint(env);

    int sorted_unit_id[MAX_PERF_UNITS];
    for (int i = 0; i < env->unit_num; ++i) {
        sorted_unit_id[i] = i;
    }

    for (int i = 0; i < env->unit_num; ++i) {
        for (int j = i + 1; j < env->unit_num; ++j) {
            if (env->total_clocks[sorted_unit_id[i]] < env->total_clocks[sorted_unit_id[j]]) {
                int tmp = sorted_unit_id[i];
                sorted_unit_id[i] = sorted_unit_id[j];
                sorted_unit_id[j] = tmp;
            }
        }
    }

    printf("%-20s|%-20s|%20s\n", "unit_name", "clocks", "percentage");
    printf("\n");
    for (int i = 0; i < env->unit_num; ++i) {
        int unit_id = sorted_unit_id[i];
        printf("%-20s|%-20lu|%20.2f%%\n", env->unit_name[unit_id], env->total_clocks[unit_id],
 env->total_clocks[unit_id] * 100.0 / env->clocks);
    }
}

#endif
