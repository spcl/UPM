#include <linux/errno.h>
#include <linux/mm.h>
#include <linux/fs.h>
#include <linux/mman.h>
#include <linux/sched.h>
#include <linux/sched/mm.h>
#include <linux/sched/coredump.h>
#include <linux/rwsem.h>
#include <linux/pagemap.h>
#include <linux/rmap.h>
#include <linux/spinlock.h>
#include <linux/slab.h>
#include <linux/memory.h>
#include <linux/hashtable.h>
#include <linux/oom.h>
#include <linux/xxhash.h>
#include <linux/compiler.h>
#include <linux/page_ref.h>
#include <linux/time64.h>

#include <asm/tlbflush.h>
#include <asm/pgtable.h>
#include <asm/spinlock.h>
#include <asm/current.h>

#include "internal.h"

#define DEFAULT_XXHASH_SEED 0
#define USM_MAX_SIZE (200 * 1024 * 1024)
#define HASH_INDEX_SIZE_COEFFICIENT 1.3

static int nr_pages_added;
static int nr_pages_replaced;

struct page_node {
    u64 hash_value;
    unsigned long addr;
	struct mm_struct *mm;
    struct page *page;
    struct hlist_node hlist_link;
};

struct rmap_node {
    u64 old_hash_value;
    unsigned long addr;
    struct mm_struct *mm;
    int pid;
    struct hlist_node hlist_link;
};

static int npage_hash;
static int nrmap_hash;
static struct hlist_head *page_hash_table;
static struct hlist_head *rmap_hash_table;
static struct kmem_cache *page_node_cache;
static struct kmem_cache *rmap_node_cache;

static DEFINE_SPINLOCK(page_hash_lock);
static DEFINE_SPINLOCK(rmap_hash_lock);

static void iterate_hash_table(void)
{
    int i;
    struct page_node *node;
    struct rmap_node *rmap_node;
    int num_page_table = 0;
    int num_rmap_table = 0;
    for (i = 0, node = NULL; node == NULL && i < npage_hash; i++)
        hlist_for_each_entry(node, &page_hash_table[i], hlist_link) {
            num_page_table++;
        }
    // pr_info("%d elements remains in page hash table\n", num_page_table);
    for (i = 0, rmap_node = NULL; rmap_node == NULL && i < nrmap_hash; i++)
        hlist_for_each_entry(rmap_node, &rmap_hash_table[i], hlist_link) {
            num_rmap_table++;
        }
    // pr_info("%d elements remains in rmap hash table\n", num_rmap_table);
}

// static void iterate_remove_hash_table(void)
// {
//     int i;
//     struct page_node *node;
//     struct hlist_node *tmp1;
//     struct rmap_node *rmap_node;
//     struct hlist_node *tmp2;
//     pr_info("usm: REMOVE page hash table and rmap hash table...\n");
//     for (i = 0, node = NULL; node == NULL && i < npage_hash; i++)
//         hlist_for_each_entry_safe(node, tmp1, &page_hash_table[i], hlist_link) {
//             hash_del(&node->hlist_link);
//         }

//     for (i = 0, rmap_node = NULL; rmap_node == NULL && i < nrmap_hash; i++)
//         hlist_for_each_entry_safe(rmap_node, tmp2, &rmap_hash_table[i], hlist_link) {
//             hash_del(&rmap_node->hlist_link);
//         }
// }

static inline int page_hash_table_init(void)
{
    npage_hash =  HASH_INDEX_SIZE_COEFFICIENT * USM_MAX_SIZE / PAGE_SIZE;
    page_hash_table = kzalloc(npage_hash * sizeof(struct hlist_head), GFP_KERNEL);
	if (!page_hash_table) {
        pr_err("usm: cannot allocate space for page hash table, quiting\n");
        return -1;
    }
    return 0;
}

static inline int rmap_hash_table_init(void)
{
    nrmap_hash =  HASH_INDEX_SIZE_COEFFICIENT * USM_MAX_SIZE / PAGE_SIZE;
    rmap_hash_table = kzalloc(nrmap_hash * sizeof(struct hlist_head), GFP_KERNEL);
	if (!rmap_hash_table) {
        pr_err("usm: cannot allocate space for rmap hash table, quiting\n");
        return -1;
    }
    return 0;
}

static inline void page_hash_table_free(void)
{
	int i;
	struct hlist_node *tmp;
    struct page_node *cur;
    for (i = 0; i < npage_hash; i++) {
        hlist_for_each_entry_safe(cur, tmp, &page_hash_table[i], hlist_link) {
            hash_del(&cur->hlist_link);
            kmem_cache_free(page_node_cache, cur);
        }
    }
    kfree(page_hash_table);
}

static inline void rmap_hash_table_free(void)
{
	int i;
	struct hlist_node *tmp;
    struct rmap_node *cur;
    for (i = 0; i < nrmap_hash; i++) {
        hlist_for_each_entry_safe(cur, tmp, &rmap_hash_table[i], hlist_link) {
            hash_del(&cur->hlist_link);
            kmem_cache_free(rmap_node_cache, cur);
        }
    }
    kfree(rmap_hash_table);
}

static inline int page_node_cache_init(void)
{
    page_node_cache = kmem_cache_create("page_node_cache", sizeof(struct page_node),
                                        sizeof(struct page_node), SLAB_HWCACHE_ALIGN, NULL);
    if (!page_node_cache) {
        pr_err("usm: cannot allocate space for page node cache, quiting\n");
        return -1;
    }
    return 0;
}

static inline int rmap_node_cache_init(void)
{
    rmap_node_cache = kmem_cache_create("rmap_node_cache", sizeof(struct rmap_node),
                                        sizeof(struct rmap_node), SLAB_HWCACHE_ALIGN, NULL);
    if (!rmap_node_cache) {
        pr_err("usm: cannot allocate space for rmap node cache, quiting\n");
        return -1;
    }
    return 0;
}

static inline void page_node_cache_free(void)
{
    kmem_cache_destroy(page_node_cache);
}

static inline void rmap_node_cache_free(void)
{
    kmem_cache_destroy(rmap_node_cache);
}

static inline struct page_node *alloc_page_node(void)
{
	void *node;
    node = kmem_cache_alloc(page_node_cache, GFP_KERNEL);
    return (struct page_node *)node;
}

static inline struct rmap_node *alloc_rmap_node(void)
{
    void *node;
    node = kmem_cache_alloc(rmap_node_cache, GFP_KERNEL);
    return (struct rmap_node *)node;
}

static inline void free_page_node(struct page_node* page_node)
{
    kmem_cache_free(page_node_cache, page_node);
}

static inline void free_rmap_node(struct rmap_node* rmap_node)
{
    kmem_cache_free(rmap_node_cache, rmap_node);
}

static inline u64 get_hash_value(struct page *page)
{
	u64 hash;
	void *addr = kmap_atomic(page);
	hash = xxh64(addr, PAGE_SIZE, DEFAULT_XXHASH_SEED);
	kunmap_atomic(addr);
	return hash;
}

static inline int get_hash_index(u64 hash_value)
{
	return hash_value % npage_hash;
}

static int add_page_to_hash_table(u64 hash_value, int hash_index, struct page *page,
                unsigned long addr, struct mm_struct *mm, u64 *spin_time)
{
    u64 sec = 0;
    u64 nsec = 0;
    struct timespec64 spin_before;
    struct timespec64 spin_after;

    struct page_node *new_page_node = alloc_page_node();
    struct rmap_node *new_rmap_node = alloc_rmap_node();
    if (!new_page_node || !new_rmap_node) {
        pr_err("usm: can't add page %px at addr %ld to hash table\n", page, addr);
        return -ENOMEM;
    }
    // pr_info("Adding page %px, pid %d, mm %px to the page tables...\n", page, current->pid, mm);
    new_page_node->hash_value = hash_value;
    new_page_node->page = page;
    new_page_node->mm = mm;
    new_page_node->addr = addr;
    ktime_get_ts64(&spin_before);
    spin_lock(&page_hash_lock);
    ktime_get_ts64(&spin_after);
    sec = spin_after.tv_sec - spin_before.tv_sec;
    nsec = sec * 1000000000 + spin_after.tv_nsec - spin_before.tv_nsec;
    *spin_time += nsec;
    hlist_add_head(&new_page_node->hlist_link, &page_hash_table[hash_index]);

    ktime_get_ts64(&spin_before);
    spin_unlock(&page_hash_lock);
    ktime_get_ts64(&spin_after);
    sec = spin_after.tv_sec - spin_before.tv_sec;
    nsec = sec * 1000000000 + spin_after.tv_nsec - spin_before.tv_nsec;
    *spin_time += nsec;

    new_rmap_node->old_hash_value = hash_value;
    new_rmap_node->addr = addr;
    new_rmap_node->mm = mm;
    new_rmap_node->pid = current->pid;
    ktime_get_ts64(&spin_before);
    spin_lock(&rmap_hash_lock);
    ktime_get_ts64(&spin_after);
    sec = spin_after.tv_sec - spin_before.tv_sec;
    nsec = sec * 1000000000 + spin_after.tv_nsec - spin_before.tv_nsec;
    *spin_time += nsec;
    hlist_add_head(&new_rmap_node->hlist_link, &rmap_hash_table[(addr + current->pid) % nrmap_hash]);
    ktime_get_ts64(&spin_before);
    spin_unlock(&rmap_hash_lock);
    ktime_get_ts64(&spin_after);
    sec = spin_after.tv_sec - spin_before.tv_sec;
    nsec = sec * 1000000000 + spin_after.tv_nsec - spin_before.tv_nsec;
    *spin_time += nsec;

    nr_pages_added++;
    return 0;
}

static int write_protect_page(struct page *page, struct vm_area_struct *vma)
{
    struct page_vma_mapped_walk pvmw = {
        .page = page,
        .vma = vma,
    };

    /* check if the pages are still valid */
    pvmw.address = page_address_in_vma(page, vma);
    if (pvmw.address == -EFAULT) {
        return -1;
    }
    if (!page_vma_mapped_walk(&pvmw)) {
        return -1;
    }

    *pvmw.pte = pte_wrprotect(*pvmw.pte);
    page_vma_mapped_walk_done(&pvmw);
    return 0;
}

/* returns 0 if identical, 1 otherwise */
static inline int do_byte_by_byte_cmp(struct page *page1, struct page *page2)
{
    int r;
    void *addr1 = kmap_atomic(page1);
    void *addr2 = kmap_atomic(page2);
    r = memcmp(addr1, addr2, PAGE_SIZE);
    kunmap_atomic(addr2);
    kunmap_atomic(addr1);
    return r;
}

static void revert_write_protect(struct page *page, struct vm_area_struct *vma)
{
    pte_t *ptep;
    spinlock_t *ptl;
    unsigned long addr;

    addr = page_address_in_vma(page, vma);
    if (vma->vm_flags & VM_WRITE) {
        ptep = get_locked_pte(vma->vm_mm, addr, &ptl);
        *ptep = pte_mkwrite(*ptep);
    }
    pte_unmap_unlock(ptep, ptl);
}

/* replace page1 pointed by addr1 with page2 */
static int replace_page(struct page *page1, struct page *page2, 
                    struct vm_area_struct *vma1, struct vm_area_struct *vma2,
                    unsigned long addr1, unsigned long addr2)
{
    pmd_t *pmd1;
    pmd_t *pmd2;
    pte_t *ptep1;
    pte_t *ptep2;
    spinlock_t *ptl1;
    spinlock_t *ptl2;
    pgprot_t pgprot;  // page protection bits of page1
    struct mm_struct *mm1 = vma1->vm_mm;
    struct mm_struct *mm2 = vma2->vm_mm;
    unsigned long address1;
    unsigned long address2;
    int err = 0;
    
    /* 
     * Check if address still points to the same struct page.
     * It's possible that user changes the content of the write-protected
     * pages before the merging, the page fault handler will allocate a
     * new physical page for the virtual address, in this case the addr
     * will match to a new page descriptor.
     */
    address1 = page_address_in_vma(page1, vma1);
    if (address1 != addr1) {
        return -EFAULT;
    }
    pmd1 = mm_find_pmd(mm1, addr1);
    ptep1 = pte_offset_map_lock(mm1, pmd1, addr1, &ptl1);

    address2 = page_address_in_vma(page2, vma2);
    if (address2 != addr2) {
        err = -EFAULT;
        goto out;
    }
    pmd2 = mm_find_pmd(mm2, addr2);
    if (mm1 == mm2 && pmd1 == pmd2) {
        ptep2 = pte_offset_map(pmd2, addr2);
    }
    else {
        ptep2 = pte_offset_map_lock(mm2, pmd2, addr2, &ptl2);
    }

    pgprot = pte_pgprot(*ptep1);

    flush_cache_page(vma1, addr1, pte_pfn(*ptep1));
    /* This includes the flushing of TLB */
    ptep_clear_flush(vma1, addr1, ptep1);
    set_pte_at(mm1, addr1, ptep1, mk_pte(page2, pgprot));

    if (PageAnon(page1)) {
        page_add_anon_rmap(page2, vma1, addr1, false);
    }
    else {
        page_add_file_rmap(page2, false);
    }

    page_remove_rmap(page1, false);

	put_page(page1);
    get_page(page2);

    if (mm1 != mm2 || pmd1 != pmd2) {
	    pte_unmap_unlock(ptep2, ptl2);
    }
out:
    pte_unmap_unlock(ptep1, ptl1);
    if (err == 0) {
        // pr_info("usm: replace page succeeds\n");
    }
    return err;
}

/* return 1 if find identical page and get them merged, 0 otherwise */
static int search_hash_table(u64 hash_value, int hash_index, struct page *page, 
                    unsigned long addr, struct mm_struct *mm, 
                    struct vm_area_struct *vma, unsigned long *vm_flags, u64 *spin_time, u64 *add_time, u64 *replace_time)
{
    u64 sec = 0;
    u64 nsec = 0;
    int r, result;
    struct page_node *cur_page_node;
    spinlock_t *ptl;
    pte_t *ptep;
    struct hlist_head *bucket;
    struct hlist_node *n;
    struct hlist_node *node;
    struct rmap_node *rmap_node;
    struct page *hash_page;
    struct vm_area_struct *cur_page_node_vma;
    struct mm_struct *cur_page_node_mm;
    struct timespec64 add_before;
    struct timespec64 add_after;
    struct timespec64 replace_before;
    struct timespec64 replace_after;
    struct timespec64 spin_before;
    struct timespec64 spin_after;
    // pr_info("search madvised page: %px at addr %ld, vma %px, mm %px in hash tables\n",\
            page, addr, vma, mm);

    /* iterate through each page in the hash table has the same hash_index */
    
    ktime_get_ts64(&spin_before);
    spin_lock(&page_hash_lock);
    ktime_get_ts64(&spin_after);
    sec = spin_after.tv_sec - spin_before.tv_sec;
    nsec = sec * 1000000000 + spin_after.tv_nsec - spin_before.tv_nsec;
    *spin_time += nsec;
    hlist_for_each_entry_safe(cur_page_node, node, &page_hash_table[hash_index], hlist_link){
        if (cur_page_node->hash_value == hash_value) {
            /* check if cur_page_node is still present */
            ptep = get_locked_pte(cur_page_node->mm, cur_page_node->addr, &ptl);
            if (!ptep || !pte_present(*ptep)) {
                pte_unmap_unlock(ptep, ptl);
                goto delete_hash_table_node;
            }
            pte_unmap_unlock(ptep, ptl);
            /* 
             * Check if cur_page_node's content has been changed
             * since it is added to the hash table.
             */
            hash_page = pfn_to_page(pte_pfn(*ptep));
            if (unlikely(hash_page == NULL)) {
                goto delete_hash_table_node;
            }
            if (get_hash_value(hash_page) != cur_page_node->hash_value) {
                goto delete_hash_table_node;
            }
            if (hash_page != cur_page_node->page) {
                goto delete_hash_table_node;
            }
            /* 
            * Check if cur_page_node is of the same page type
            * (both anonymous or both file-backed) as the madvised page.
            */
            if (PageAnon(hash_page) ^ PageAnon(page)) {
                continue;
            }
            /* check if current checking address already points to the same physical page */
            if (unlikely(page_to_pfn(page) == page_to_pfn(hash_page))) {
                spin_unlock(&page_hash_lock);
                return 0;
            }

            cur_page_node_mm = cur_page_node->mm;
            cur_page_node_vma = find_vma(cur_page_node_mm, cur_page_node->addr);
            // down_read(&mm->mmap_sem);
            // down_read(&cur_page_node_mm->mmap_sem);
            /* lock the pages in order to mlock them */
            lock_page(page);
            lock_page(hash_page);
            mlock_vma_page(page);
            mlock_vma_page(hash_page);

            /*
             * Write-protect both pages.
             * If write-protect fails, continue to the next page in the hash table.
             */
            r = write_protect_page(page, vma);
            if (r < 0) {
                // pr_info("usm: wrprotect madvised page failed\n");
                goto r1;
            }
            r = write_protect_page(hash_page, cur_page_node_vma);
            if (r < 0) {
                // pr_info("usm: wrprotect hash_page failed\n");
                goto r2;
            }
            if (do_byte_by_byte_cmp(page, cur_page_node->page) == 0) {
                // pr_info("found page %px fwith the same content with page %px\n",\
                        page, cur_page_node->page);
                /* merge the pages */
                ktime_get_ts64(&replace_before);
                result = replace_page(page, cur_page_node->page, vma, 
                                cur_page_node_vma, addr, cur_page_node->addr);
                ktime_get_ts64(&replace_after);
                sec = replace_after.tv_sec - replace_before.tv_sec;
                nsec = sec * 1000000000 + replace_after.tv_nsec - replace_before.tv_nsec;
                *replace_time += nsec;
                if (result < 0) {
                    /*
                     * If replace_page fails, there won't be another page in the hash
                     * tablefre that has the same content with the madvised page (if there is,
                     * it should be merged with hash_page), so there's
                     * no need continuing searching the hash table, we can directly go to
                     * the next madvised page.
                     */
                    // pr_info("usm: replace page %px by page %px failed\n", page, cur_page_node->page);
                    r = -2;
                    goto r2;
                }
                else {
                    /* replace succeeds */
                    r = 0;
                    // pr_info("madvised page %px mapcount %d, refcount %d\n", page, page_mapcount(page), page_ref_count(page));
                    // pr_info("hash page %px, mapcount %d, refcount %d\n", hash_page, page_mapcount(page), page_ref_count(page));
                    nr_pages_replaced++;
                    goto unlock;
                }
            }
            /* byte-by-byte comparison fails, go to the next page in the hash table */
            else {
                r = -1;
                goto r2;
            }
        r2:
            revert_write_protect(hash_page, cur_page_node_vma);
        r1:
            revert_write_protect(page, vma);
            goto unlock;

        delete_hash_table_node:
            bucket = &rmap_hash_table[(cur_page_node->addr + current->pid) % nrmap_hash];
            spin_lock(&rmap_hash_lock);
            hlist_for_each_entry_safe(rmap_node, n, bucket, hlist_link) {
                if (rmap_node->addr == cur_page_node->addr &&\
                    rmap_node->mm == cur_page_node->mm) {
                    // spin_lock(&rmap_hash_lock);
                    hash_del(&rmap_node->hlist_link);
                    // spin_unlock(&rmap_hash_lock);
                    free_rmap_node(rmap_node);
                    break;
                }
            }
            spin_unlock(&rmap_hash_lock);
            hash_del(&cur_page_node->hlist_link);
            free_page_node(cur_page_node);
            continue;

        unlock:
            munlock_vma_page(hash_page);
            munlock_vma_page(page);
            unlock_page(hash_page);
            unlock_page(page);
            // up_read(&cur_page_node_mm->mmap_sem);
            // up_read(&mm->mmap_sem);
            if (r == -1) {
                continue;
            }
            else {
                ktime_get_ts64(&spin_before);
                spin_unlock(&page_hash_lock);
                ktime_get_ts64(&spin_after);
                sec = spin_after.tv_sec - spin_before.tv_sec;
                nsec = sec * 1000000000 + spin_after.tv_nsec - spin_before.tv_nsec;
                *spin_time += nsec;
                return r;
            }
        }
    }
    ktime_get_ts64(&spin_before);
    spin_unlock(&page_hash_lock);
    ktime_get_ts64(&spin_after);
    sec = spin_after.tv_sec - spin_before.tv_sec;
    nsec = sec * 1000000000 + spin_after.tv_nsec - spin_before.tv_nsec;
    *spin_time += nsec;

    /*
     * If the function haven't returned yet, this means no matching page is found,
     * add this page to the hash table.
     */
    ktime_get_ts64(&add_before);
    r = add_page_to_hash_table(hash_value, hash_index, page, addr, mm, spin_time);
    ktime_get_ts64(&add_after);
    sec = add_after.tv_sec - add_before.tv_sec;
    nsec = sec * 1000000000 + add_after.tv_nsec - add_before.tv_nsec;
    *add_time += nsec;
    
    if (r == 0) {
        current->flags |= PF_USM;
    }
    return r;
}

/*	Now we assume that there is only advise behavior, no un-adivse one.
 *	
 */
int usm_madvise(struct vm_area_struct *vma, unsigned long start,
		unsigned long end, int advice, unsigned long *vm_flags)
{
	unsigned long cur_addr;
    struct page *cur_page;
    struct mm_struct *mm;
    struct rmap_node *old_rmap_node = NULL;
    struct page_node *old_page_node = NULL;
    struct hlist_head *rmap_bucket;
    struct hlist_head *page_bucket;
    struct hlist_node *cur_node, *n;
    u64 hash_value;
    int hash_index;
    u64 search_time = 0;
    u64 follow_page_time = 0;
    u64 hash_time = 0;
    u64 madvise_time = 0;
    u64 check_time = 0;
    u64 spin_time = 0;
    u64 add_time = 0;
    u64 replace_time = 0;
    u64 sec = 0;
    u64 nsec = 0;
    struct timespec64 madvise_before;
    struct timespec64 madvise_after;
    struct timespec64 hash_before;
    struct timespec64 hash_after;
    struct timespec64 follow_before;
    struct timespec64 follow_after;
    struct timespec64 search_before;
    struct timespec64 search_after;
    struct timespec64 check_before;
    struct timespec64 check_after;
    struct timespec64 spin_before;
    struct timespec64 spin_after;

    // int counter = 0;
    ktime_get_ts64(&madvise_before);

    if (*vm_flags & (VM_SHARED  | VM_MAYSHARE   |
				 VM_PFNMAP    | VM_IO      | VM_DONTEXPAND |
				 VM_HUGETLB | VM_MIXEDMAP))
		return 0;

    mm = vma->vm_mm;
    if (mm) {
        mmgrab(mm);
    }

    // pr_info("Before going into usm, iterating and remove page tables:\n");
	for (cur_addr = start; cur_addr < end; cur_addr += PAGE_SIZE) {
        // down_read(&mm->mmap_sem);
        // pr_info("counter: %d\n", ++counter);
        // pr_info("checking cur_addr %ld\n", cur_addr);
        /* get the struct page at address cur_addr */
        ktime_get_ts64(&follow_before);
        cur_page = follow_page(vma, cur_addr, FOLL_GET);
        ktime_get_ts64(&follow_after);
        sec = follow_after.tv_sec - follow_before.tv_sec;
        nsec = sec * 1000000000 + follow_after.tv_nsec - follow_before.tv_nsec;
        follow_page_time += nsec;
        // up_read(&mm->mmap_sem);
        if (cur_page == NULL) {
            // pr_info("usm: no present page at address %lx\n", cur_addr);
            continue;
        }
        // pr_info("cur_page: %px, mapcount %d, refcount %d\n", cur_page, page_mapcount(cur_page), page_ref_count(cur_page));

        ktime_get_ts64(&hash_before);
        hash_value = get_hash_value(cur_page);
        ktime_get_ts64(&hash_after);
        sec = hash_after.tv_sec - hash_before.tv_sec;
        nsec = sec * 1000000000 + hash_after.tv_nsec - hash_before.tv_nsec;
        hash_time += nsec;
        hash_index = get_hash_index(hash_value);

        /* check if the page is already stored by USM */
        rmap_bucket = &(rmap_hash_table[(cur_addr + current->pid) % nrmap_hash]);
        ktime_get_ts64(&check_before);
        ktime_get_ts64(&spin_before);
        spin_lock(&rmap_hash_lock);
        ktime_get_ts64(&spin_after);
        sec = spin_after.tv_sec - spin_before.tv_sec;
        nsec = sec * 1000000000 + spin_after.tv_nsec - spin_before.tv_nsec;
        spin_time += nsec;

        hlist_for_each_entry_safe(old_rmap_node, n, rmap_bucket, hlist_link) {
            if (old_rmap_node->addr == cur_addr && old_rmap_node->mm == mm) {
                if (old_rmap_node->old_hash_value != hash_value) {
                    hash_del(&old_rmap_node->hlist_link);
                    page_bucket = &(page_hash_table[old_rmap_node->old_hash_value % npage_hash]);
                    spin_unlock(&rmap_hash_lock);
                    spin_lock(&page_hash_lock);
                    hlist_for_each_entry_safe(old_page_node, cur_node, page_bucket, hlist_link) {
                        if (old_page_node->hash_value == old_rmap_node->old_hash_value) {
                            hash_del(&old_page_node->hlist_link);
                            break;
                        }
                    }
                    spin_unlock(&page_hash_lock);
                    spin_lock(&rmap_hash_lock);
                    free_page_node(old_page_node);
                    free_rmap_node(old_rmap_node);
                    goto search_hash_table;
                }
                else {
                    /* the page is already added to USM, no need to do anything */
                    spin_unlock(&rmap_hash_lock);
                    ktime_get_ts64(&check_after);
                    sec = check_after.tv_sec - check_before.tv_sec;
                    nsec = sec * 1000000000 + check_after.tv_nsec - check_before.tv_nsec;
                    check_time += nsec;
                    goto next;
                }
            }
        }
    search_hash_table:
        ktime_get_ts64(&spin_before);
        spin_unlock(&rmap_hash_lock);
        ktime_get_ts64(&spin_after);
        sec = spin_after.tv_sec - spin_before.tv_sec;
        nsec = sec * 1000000000 + spin_after.tv_nsec - spin_before.tv_nsec;
        spin_time += nsec;

        ktime_get_ts64(&check_after);
        sec = check_after.tv_sec - check_before.tv_sec;
        nsec = sec * 1000000000 + check_after.tv_nsec - check_before.tv_nsec;
        check_time += nsec;
        // counter += 1;
        ktime_get_ts64(&search_before);
        search_hash_table(hash_value, hash_index, cur_page, cur_addr, mm, vma, vm_flags, &spin_time, &add_time, &replace_time);
        ktime_get_ts64(&search_after);
        sec = search_after.tv_sec - search_before.tv_sec;
        nsec = sec * 1000000000 + search_after.tv_nsec - search_before.tv_nsec;
        search_time += nsec;
    next:
        // iterate_hash_table();
        put_page(cur_page);
        continue;
    }

    mmdrop(mm);
    ktime_get_ts64(&madvise_after);
    sec = madvise_after.tv_sec - madvise_before.tv_sec;
    nsec = sec * 1000000000 + madvise_after.tv_nsec - madvise_before.tv_nsec;
    madvise_time += nsec;
    // pr_info("task pid %d has madvised address %lx, %ld pages, added %d pages to USM, replaced %d pages\n",\
            current->pid, start, (end-start)/4096, nr_pages_added, nr_pages_replaced);
    nr_pages_added = 0;
    nr_pages_replaced = 0;
    pr_info("search_time %lld add_time %lld replace_time %lld spin_time %lld, follow_time %lld hash_time %lld check_time %lld  madvise_time %lld ",
            search_time / 1000, add_time / 1000, replace_time / 1000, spin_time / 1000, follow_page_time / 1000, hash_time / 1000, check_time / 1000, madvise_time / 1000);
    // pr_info("enter into search_hash_time %d times\n", counter);
    return 0;
}

/* search in the mm the pages stored in the hash tables, and delete then */
void usm_exit(int pid, struct mm_struct *mm)
{
    struct page_node *cur_page_node;
    struct hlist_node *tmp_page_node;
    struct rmap_node *cur_rmap_node;
    struct hlist_node *tmp_rmap_node;
    u64 hash_val;
    int i;
    u64 exit_time = 0;
    u64 sec = 0;
    u64 nsec = 0;
    struct timespec64 exit_before;
    struct timespec64 exit_after;

    ktime_get_ts64(&exit_before);

    // pr_info("usm_exit for mm %px, pid %d\n", mm, pid);
    for (i = 0, cur_rmap_node = NULL; i < nrmap_hash; i++) {
        hlist_for_each_entry_safe(cur_rmap_node, tmp_rmap_node, &rmap_hash_table[i], hlist_link) {
            if (cur_rmap_node->pid != pid)
                continue;
            hash_val = cur_rmap_node->old_hash_value;
            // pr_info("find rmap page addr %ld, mm %px, pid %d\n", cur_rmap_node->addr, cur_rmap_node->mm, pid);
            hlist_for_each_entry_safe(cur_page_node, tmp_page_node,\
                &page_hash_table[hash_val % npage_hash], hlist_link) {
                if (cur_page_node->hash_value != hash_val ||\
                    cur_page_node->mm != cur_rmap_node->mm || cur_page_node->addr != cur_rmap_node->addr)
                    continue;
                // pr_info("find page %px, addr %ld, deleting it...\n", cur_page_node->page, cur_page_node->addr);
                hash_del(&cur_page_node->hlist_link);
                free_page_node(cur_page_node);
                break;
            }
            hash_del(&cur_rmap_node->hlist_link);
            free_rmap_node(cur_rmap_node);
            continue;
        }
    }
    current->flags &= ~PF_USM;

    ktime_get_ts64(&exit_after);
    sec = exit_after.tv_sec - exit_before.tv_sec;
    nsec = sec * 1000000000 + exit_after.tv_nsec - exit_before.tv_nsec;
    exit_time += nsec;
    pr_info("usm_exit time: %lld\n", exit_time);
    // iterate_hash_table();
}

static int __init usm_init(void)
{
    int error = 0;
    pr_info("Entering usm_init...\n");
    /* initialize the data structures */
    if (page_hash_table_init() < 0) {
        pr_err("usm: page_hash_table_init failed\n");
        return -ENOMEM;
    }
    if (page_node_cache_init() < 0) {
        error = -ENOMEM;
        pr_err("usm: page_node_cache_init failed\n");
        goto quit_free1;
    }
    if (rmap_hash_table_init() < 0) {
        error = -ENOMEM;
        pr_err("usm: rmap_hash_table_init failed\n");
        goto quit_free2;
    }
    if (rmap_node_cache_init() < 0) {
        error = -ENOMEM;
        pr_err("usm: rmap_node_cache_init failed\n");
        goto quit_free3;
    }

    return 0;

quit_free3:
    rmap_hash_table_free();
quit_free2:
    page_node_cache_free();
quit_free1:
    page_hash_table_free();

    if (error < 0) {
        pr_err("usm init failed\n");
    }
    return error;
}
subsys_initcall(usm_init);
