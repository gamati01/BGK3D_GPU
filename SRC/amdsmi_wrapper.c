#ifdef ENERGY
/*
 * AMD SMI counterpart of nvml_wrapper.c
 *
 * Provides the very same C symbols used by the Fortran side
 * (nvml_interface_module.F90 / bgk3d.F90), so the Fortran code does not
 * need any change: only the wrapper object that gets linked changes.
 *
 *   get_gpu_power_index(idx)        -> instantaneous GPU power  [milliWatt]
 *   get_gpu_energy_mJ_u64(idx,*e)   -> cumulative GPU energy    [milliJoule]
 *
 * NVML equivalences:
 *   nvmlDeviceGetPowerUsage              <-> amdsmi_get_power_info
 *   nvmlDeviceGetTotalEnergyConsumption  <-> amdsmi_get_energy_count
 *
 * Build (see Makefile rule), e.g.:
 *   clang -O3 -I/opt/rocm/include -c amdsmi_wrapper.c
 *   link with: -L/opt/rocm/lib -lamd_smi
 */
#include <stdint.h>
#include <amd_smi/amdsmi.h>

/* Max number of GPUs we are willing to enumerate */
#define MAX_GPUS 64

/* Flat list of GPU processor handles, indexed the same way NVML indexes
 * devices (0,1,2,...). Built once, lazily. */
static amdsmi_processor_handle g_dev[MAX_GPUS];
static uint32_t                g_ndev   = 0;
static int                     g_inited = 0;

/* Initialize AMD SMI and build the flat device list.
 * Returns 0 on success, negative on failure. */
static int amdsmi_setup(void) {
    if (g_inited) return 0;

    if (amdsmi_init(AMDSMI_INIT_AMD_GPUS) != AMDSMI_STATUS_SUCCESS) return -1;

    uint32_t socket_count = 0;
    if (amdsmi_get_socket_handles(&socket_count, NULL) != AMDSMI_STATUS_SUCCESS) return -1;
    if (socket_count == 0) return -1;

    amdsmi_socket_handle sockets[MAX_GPUS];
    if (socket_count > MAX_GPUS) socket_count = MAX_GPUS;
    if (amdsmi_get_socket_handles(&socket_count, sockets) != AMDSMI_STATUS_SUCCESS) return -1;

    g_ndev = 0;
    for (uint32_t s = 0; s < socket_count && g_ndev < MAX_GPUS; ++s) {
        uint32_t pcount = 0;
        if (amdsmi_get_processor_handles(sockets[s], &pcount, NULL) != AMDSMI_STATUS_SUCCESS)
            continue;
        if (pcount == 0) continue;

        amdsmi_processor_handle procs[MAX_GPUS];
        if (pcount > MAX_GPUS) pcount = MAX_GPUS;
        if (amdsmi_get_processor_handles(sockets[s], &pcount, procs) != AMDSMI_STATUS_SUCCESS)
            continue;

        for (uint32_t p = 0; p < pcount && g_ndev < MAX_GPUS; ++p)
            g_dev[g_ndev++] = procs[p];
    }

    if (g_ndev == 0) return -1;

    g_inited = 1;
    return 0;
}

/* Instantaneous power of GPU 'idx' in milliWatt (to mirror NVML, which
 * reports milliWatt). AMD SMI reports socket power in Watt. Returns a
 * negative value on error. */
int get_gpu_power_index(int idx) {
    if (amdsmi_setup() != 0) return -1;
    if (idx < 0 || (uint32_t)idx >= g_ndev) return -2;

    amdsmi_power_info_t info;
    if (amdsmi_get_power_info(g_dev[idx], &info) != AMDSMI_STATUS_SUCCESS) return -3;

    /* Prefer the instantaneous reading; fall back to the average for the
     * cards that only populate average_socket_power (Navi / MI200 and older).
     * Unsupported members are reported as UINT32_MAX by AMD SMI. */
    uint32_t p_w = info.current_socket_power;
    if (p_w == 0u || p_w == 0xFFFFFFFFu) p_w = info.average_socket_power;
    if (p_w == 0xFFFFFFFFu) return -4;

    return (int)(p_w * 1000u);   /* W -> mW */
}

/* Cumulative energy of GPU 'idx' in milliJoule (to mirror NVML's
 * nvmlDeviceGetTotalEnergyConsumption). AMD SMI returns a raw accumulator
 * plus a resolution; energy[uJ] = accumulator * counter_resolution.
 * Returns 0 on success, negative on error. */
int get_gpu_energy_mJ_u64(int idx, unsigned long long *e_mJ) {
    if (e_mJ == 0) return -5;
    if (amdsmi_setup() != 0) return -1;
    if (idx < 0 || (uint32_t)idx >= g_ndev) return -2;

    uint64_t accumulator      = 0;
    float    counter_res      = 0.0f;
    uint64_t timestamp        = 0;
    if (amdsmi_get_energy_count(g_dev[idx], &accumulator, &counter_res, &timestamp)
            != AMDSMI_STATUS_SUCCESS)
        return -3;

    /* accumulator * counter_resolution -> microJoule, then uJ -> mJ */
    double energy_uJ = (double)accumulator * (double)counter_res;
    *e_mJ = (unsigned long long)(energy_uJ / 1000.0);

    return 0;
}
#endif
