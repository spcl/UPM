#ifndef __LINUX_USM_H
#define __LINUX_USM_H
/*
 * Memory merging support.
 *
 * This code enables dynamic sharing of identical pages found in different
 * memory areas, even if they are not shared by fork().
 */

// #include <linux/bitops.h>
// #include <linux/mm.h>
// #include <linux/pagemap.h>
// #include <linux/rmap.h>
// #include <linux/sched.h>
// #include <linux/sched/coredump.h>


int usm_madvise(struct vm_area_struct *vma, unsigned long start,
		unsigned long end, int advice, unsigned long *vm_flags);

#endif /* __LINUX_KSM_H */
