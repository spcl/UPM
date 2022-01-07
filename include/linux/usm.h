#ifndef __LINUX_USM_H
#define __LINUX_USM_H

int usm_madvise(struct vm_area_struct *vma, unsigned long start,
		unsigned long end, int advice, unsigned long *vm_flags);

// void usm_exit(struct mm_struct *mm);
void usm_exit(int pid, struct mm_struct *mm);

#endif /* __LINUX_USM_H */
