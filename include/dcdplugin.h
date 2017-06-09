#ifndef DCDPLUGIN_H
#define DCDPLUGIN_H
#include "molfile_plugin.h"
#include "fastio.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  fio_fd fd;
  int natoms;
  int nsets;
  int setsread;
  int istart;
  int nsavc;
  double delta;
  int nfixed;
  float *x, *y, *z;
  int *freeind;
  float *fixedcoords;
  int reverse;
  int charmm;  
  int first;
  int with_unitcell;
} dcdhandle;


extern void *open_dcd_read(const char *path, const char *filetype, int *natoms);
extern int read_next_timestep(void *v, int natoms, molfile_timestep_t *ts);
extern void close_file_read(void *v);
extern int get_nframes(void *v);

#ifdef __cplusplus
}
#endif

#endif
