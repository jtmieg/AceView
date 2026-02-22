/*
 * sa.seeds.c

 * This module is part of the sortalign package
 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 * Created April 18, 2025

 * This is public.


 * This module analyses the hardware
 * to bing the aligner to the least buzy core
*/

#include <dirent.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/sysinfo.h>  // for get_nprocs_conf() alternative if needed
#include <unistd.h>       // usleep, access
#include <time.h>         // time(NULL)
#include <stdlib.h>       // random, srandom
#include <string.h>       // random, srandom

#define MAX_CPUS 256  // Adjust if your machine has more

typedef struct {
    unsigned long long user, nice, sys, idle, iowait, irq, softirq, steal, guest, guest_nice ;
} cpu_times_t ;

// Function to read /proc/stat into per-CPU times array (returns num_cpus read)

static int read_proc_stat(cpu_times_t *times)
{
    FILE *f = fopen("/proc/stat", "r") ;
    if (!f) return -1 ;

    char buf[2048] ;
    int num_cpus = 0 ;
    int cpu_id = 0 ;
    
    // Parse cpuN + the 10 counters (user nice system idle iowait irq softirq steal guest guest_nice)
    if (sscanf(buf, "cpu%d %llu %llu %llu %llu %llu %llu %llu %llu %llu %llu",
	       &cpu_id,
	       &times[cpu_id].user,
	       &times[cpu_id].nice,
	       &times[cpu_id].sys,
	       &times[cpu_id].idle,
	       &times[cpu_id].iowait,
	       &times[cpu_id].irq,
	       &times[cpu_id].softirq,
	       &times[cpu_id].steal,
	       &times[cpu_id].guest,
	       &times[cpu_id].guest_nice) == 11)
      {
	    
	if (cpu_id >= 0 && cpu_id < MAX_CPUS) {
	  if (cpu_id + 1 > num_cpus) num_cpus = cpu_id + 1 ;
	}
      }
    
  fclose(f) ;
  return num_cpus ;
}
  
/* node with least average CPU utilization */
int saBestNumactlNode (void)
{
  int best_node = 0 ;
  double min_avg_usage = -1.0 ;

  // First sample
  cpu_times_t before[MAX_CPUS] = {0} ;
  int num_cpus = read_proc_stat(before) ;
  if (num_cpus <= 1) // no choice
    return 0 ;      
  
  usleep(500000) ;  // 0.5s delay for delta (tune if needed)
  srandom(time(NULL) ^ getpid()) ;
  best_node = random() % num_cpus ; /* random fallback */

  
  // Second sample
  cpu_times_t after[MAX_CPUS] = {0} ;
  if (read_proc_stat(after) != num_cpus)
    return best_node ;
  
  
  // Compute per-CPU utilization %
  double usage[MAX_CPUS] = {0.0} ;
  for (int cpu = 0 ; cpu < num_cpus ; cpu++) {
    unsigned long long delta_user = after[cpu].user - before[cpu].user ;
    unsigned long long delta_nice = after[cpu].nice - before[cpu].nice ;
    unsigned long long delta_sys = after[cpu].sys - before[cpu].sys ;
    unsigned long long delta_idle = after[cpu].idle - before[cpu].idle ;
    unsigned long long delta_iowait = after[cpu].iowait - before[cpu].iowait ;
    unsigned long long delta_irq = after[cpu].irq - before[cpu].irq ;
    unsigned long long delta_softirq = after[cpu].softirq - before[cpu].softirq ;
    unsigned long long delta_steal = after[cpu].steal - before[cpu].steal ;

    unsigned long long delta_total = delta_user + delta_nice + delta_sys + delta_idle +
      delta_iowait + delta_irq + delta_softirq + delta_steal ;
    if (delta_total == 0) continue ;  // Avoid div0
    
    unsigned long long delta_active = delta_user + delta_nice + delta_sys +
      delta_irq + delta_softirq + delta_steal ;
    usage[cpu] = 100.0 * ((double)delta_active / (double)delta_total) ;
  }
  
  // Now loop over nodes
  DIR *d = opendir("/sys/devices/system/node") ;
  if (d) {
    struct dirent *e ;
    while ((e = readdir(d))) {
      int node ;
      if (sscanf(e->d_name, "node%d", &node) != 1) continue ;
      
      // Get cpumap for this node (your original way)
      char path[64] ;
      snprintf(path, sizeof(path), "/sys/devices/system/node/node%d/cpumap", node) ;
      FILE *f = fopen(path, "r") ;
      if (!f) continue ;
      unsigned long long map = 0 ;
      int n = fscanf(f, "%llx", &map) ;
      fclose(f) ;
      if (n != 1) continue ;
      
      // Compute avg usage for CPUs in this node
      double node_usage_sum = 0.0 ;
      int node_cpu_count = 0 ;
      for (int cpu = 0 ; cpu < num_cpus && cpu < 64 ; cpu++) {  // Assumes map fits ull
	if (map & (1ULL << cpu)) {
	  node_usage_sum += usage[cpu] ;
	  node_cpu_count++ ;
	}
      }
      if (node_cpu_count == 0) continue ;
      
      double avg_usage = node_usage_sum / node_cpu_count ;
      
      if (min_avg_usage == -1.0 || avg_usage < min_avg_usage)
	{
	  min_avg_usage = avg_usage ;
	  best_node = node ;
	}
    }
    closedir(d) ;
  } 
  
  return best_node ;
}

#ifdef JUNK

foreach ii (1 2)
  \rm -rf titi$ii ; /usr/bin/time -f "E %E U %U M %M P %P" ~/ace/bin.LINUX_4_OPT/sortalign -x Aligners/011_SortAlignG5R5/IDX.GRCh38.18.31 -i Fasta/iRefSeq38/iRefSeq38.fasta.gz --align -o titi$ii --step 5 >& titi$ii.err &
sleep 4
end

#endif
  


