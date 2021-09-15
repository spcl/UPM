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

#include <asm/tlbflush.h>
#include <asm/pgtable.h>
#include <asm/spinlock.h>

#include "internal.h"

#define DEFAULT_XXHASH_SEED 0
#define USM_MAX_SIZE (200 * 1024 * 1024)
#define HASH_INDEX_SIZE_COEFFICIENT 1.3

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
    // struct rmap_node *rmap_node;
    pr_info("usm: list all elements in the page hash table:\n");
    for (i = 0, node = NULL; node == NULL && i < npage_hash; i++)
        hlist_for_each_entry(node, &page_hash_table[i], hlist_link) {
            pr_info("page: %px, hash value: %lld, hash_index: %d, addr: %ld\n",\
                    node->page, node->hash_value, i, node->addr);
        }
/*
    pr_info("usm: list all elements in the rmap hash table:\n");
    for (i = 0, rmap_node = NULL; rmap_node == NULL && i < nrmap_hash; i++)
        hlist_for_each_entry(rmap_node, &rmap_hash_table[i], hlist_link) {
            pr_info("addr: %ld, old hash value: %lld\n",\
                    rmap_node->addr, rmap_node->old_hash_value);
        }
        */
}

static void iterate_remove_hash_table(void)
{
    int i;
    struct page_node *node;
    struct hlist_node *tmp1;
    struct rmap_node *rmap_node;
    struct hlist_node *tmp2;
    pr_info("usm: list and REMOVE page hash table:\n");
    for (i = 0, node = NULL; node == NULL && i < npage_hash; i++)
        hlist_for_each_entry_safe(node, tmp1, &page_hash_table[i], hlist_link) {
            pr_info("page: %px, hash value: %lld, hash_index: %d, addr: %ld\n",\
                    node->page, node->hash_value, i, node->addr);
            hash_del(&node->hlist_link);
        }

    pr_info("usm: list and REMOVE rmap hash table:\n");
    for (i = 0, rmap_node = NULL; rmap_node == NULL && i < nrmap_hash; i++)
        hlist_for_each_entry_safe(rmap_node, tmp2, &rmap_hash_table[i], hlist_link) {
            pr_info("addr: %ld, old hash value: %lld\n",\
                    rmap_node->addr, rmap_node->old_hash_value);
            hash_del(&rmap_node->hlist_link);
        }
}

static inline int page_hash_table_init(void)
{
    npage_hash =  HASH_INDEX_SIZE_COEFFICIENT * USM_MAX_SIZE / PAGE_SIZE;
	pr_info("usm: npage_hash: %d\n", npage_hash);
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
	pr_info("usm: nrmap_hash: %d\n", nrmap_hash);
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

static inline int add_page_to_hash_table(u64 hash_value, int hash_index, struct page *page,
                unsigned long addr, struct mm_struct *mm)
{
    struct page_node *new_page_node = alloc_page_node();
    struct rmap_node *new_rmap_node = alloc_rmap_node();
    if (!new_page_node || !new_rmap_node) {
        pr_err("usm: can't add page %px at addr %ld to hash table\n", page, addr);
        return -ENOMEM;
    }
    pr_info("Adding page at addr %ld, page %px to the page tables...\n", addr, page);
    spin_lock(&page_hash_lock);
    new_page_node->hash_value = hash_value;
    new_page_node->page = page;
    new_page_node->mm = mm;
    new_page_node->addr = addr;
    hlist_add_head(&new_page_node->hlist_link, &page_hash_table[hash_index]);
    spin_unlock(&page_hash_lock);

    spin_lock(&rmap_hash_lock);
    new_rmap_node->old_hash_value = hash_value;
    new_rmap_node->addr = addr;
    new_rmap_node->mm = mm;
    hlist_add_head(&new_rmap_node->hlist_link, &rmap_hash_table[addr % nrmap_hash]);
    spin_unlock(&rmap_hash_lock);

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
static inline int replace_page(struct page *page1, struct page *page2, 
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

    pr_info("before page_remove_rmap, page1 %px, mapcount %d, refcount %d\n",\
            page1, page_mapcount(page1), page_ref_count(page1));
    page_remove_rmap(page1, false);
    pr_info("after page_remove_rmap, page1 %px, mapcount %d, refcount %d\n",\
            page1, page_mapcount(page1), page_ref_count(page1));

    /* ksm 放在swap的也能share， 我把page mlock再memory里了*/
	// if (!page_mapped(page1))
	//     try_to_free_swap(page1);
    pr_info("before put/get page, page1 %px, mapcount %d, refcount %d\n",\
            page1, page_mapcount(page1), page_ref_count(page1));
    pr_info("before put/get page, page2 %px, mapcount %d, refcount %d\n",\
            page2, page_mapcount(page2), page_ref_count(page2));
	put_page(page1);
    get_page(page2);
    pr_info("After put/get page, page1 %px, mapcount %d, refcount %d\n",\
            page1, page_mapcount(page1), page_ref_count(page1));
    pr_info("After put/get page, page2 %px, mapcount %d, refcount %d\n",\
            page2, page_mapcount(page2), page_ref_count(page2));

    if (mm1 != mm2 || pmd1 != pmd2) {
	    pte_unmap_unlock(ptep2, ptl2);
    }
out:
    pte_unmap_unlock(ptep1, ptl1);
    if (err == 0) {
        pr_info("usm: replace page succeeds\n");
    }
    return err;
}

/* return 1 if find identical page and get them merged, 0 otherwise */
static int search_hash_table(u64 hash_value, int hash_index, struct page *page, 
                    unsigned long addr, struct mm_struct *mm, 
                    struct vm_area_struct *vma, unsigned long *vm_flags)
{
    int r;
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
    pr_info("search madvised page: %px at addr %ld, vma %px, mm %px in hash tables\n",\
            page, addr, vma, mm);

    /* iterate through each page in the hash table has the same hash_index */
    spin_lock(&page_hash_lock);
    hlist_for_each_entry_safe(cur_page_node, node, &page_hash_table[hash_index], hlist_link){
        pr_info("page in hash table: %px\n", cur_page_node);
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
            pr_info("lock page %px, is locked? %d\n", page, PageLocked(page));
            lock_page(hash_page);
            pr_info("lock hash page %px, is locked? %d\n", hash_page, PageLocked(hash_page));
            mlock_vma_page(page);
            pr_info("mlock page %px, is locked? %d\n", page, PageLocked(page));
            mlock_vma_page(hash_page);
            pr_info("mlock hash_page %px, is locked? %d\n", hash_page, PageLocked(hash_page));

            /*
             * Write-protect both pages.
             * If write-protect fails, continue to the next page in the hash table.
             */
            r = write_protect_page(page, vma);
            if (r < 0) {
                pr_info("usm: wrprotect madvised page failed\n");
                goto r1;
            }
            r = write_protect_page(hash_page, cur_page_node_vma);
            if (r < 0) {
                pr_info("usm: wrprotect hash_page failed\n");
                goto r2;
            }
            if (do_byte_by_byte_cmp(page, cur_page_node->page) == 0) {
                pr_info("found page %px fwith the same content with page %px\n",\
                        page, cur_page_node->page);
                /* merge the pages */
                if (replace_page(page, cur_page_node->page, vma, 
                                cur_page_node_vma, addr, cur_page_node->addr) < 0) {
                    /*
                     * If replace_page fails, there won't be another page in the hash
                     * tablefre that has the same content with the madvised page (if there is,
                     * it should be merged with hash_page), so there's
                     * no need continuing searching the hash table, we can directly go to
                     * the next madvised page.
                     */
                    pr_info("usm: replace page %px by page %px failed\n", page, cur_page_node->page);
                    r = -2;
                    goto r2;
                }
                else {
                    r = 0;
                    pr_info("point 5:\n");
                    pr_info("madvised page %px mapcount %d, refcount %d\n", page, page_mapcount(page), page_ref_count(page));
                    pr_info("hash page %px, mapcount %d, refcount %d\n", hash_page, page_mapcount(page), page_ref_count(page));
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
            hash_del(&cur_page_node->hlist_link);
            spin_lock(&rmap_hash_lock);
            bucket = &rmap_hash_table[cur_page_node->addr % nrmap_hash];
            hlist_for_each_entry_safe(rmap_node, n, bucket, hlist_link) {
                if (rmap_node->addr == cur_page_node->addr &&\
                    rmap_node->mm == cur_page_node->mm) {
                    hash_del(&rmap_node->hlist_link);
                    break;
                }
            }
            spin_unlock(&rmap_hash_lock);
            continue;

        unlock:
            munlock_vma_page(hash_page);
            pr_info("munlock hash page %px, is locked? %d\n", hash_page, PageLocked(hash_page));
            munlock_vma_page(page);
            pr_info("munlock page %px\n, is locked? %d", page, PageLocked(page));
            unlock_page(hash_page);
            pr_info("unlock hash page %px, is locked? %d\n", hash_page, PageLocked(hash_page));
            unlock_page(page);
            pr_info("unlock page %px, is locked? %d\n", page, PageLocked(page));
            // up_read(&cur_page_node_mm->mmap_sem);
            // up_read(&mm->mmap_sem);
            if (r == -1) {
                continue;
            }
            else {
                spin_unlock(&page_hash_lock);
                return r;
            }
        }
    }
    spin_unlock(&page_hash_lock);

    /*
     * If the function haven't returned yet, this means no matching page is found,
     * add this page to the hash table.
     */
    r = add_page_to_hash_table(hash_value, hash_index, page, addr, mm);
    if (r == 0) {
        set_bit(MMF_VM_SHAREABLE, &mm->flags);
        *vm_flags |= VM_SHAREABLE;
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
    int counter = 0;

    if (*vm_flags & (VM_SHARED  | VM_MAYSHARE   |
				 VM_PFNMAP    | VM_IO      | VM_DONTEXPAND |
				 VM_HUGETLB | VM_MIXEDMAP))
		return 0;

    mm = vma->vm_mm;
    if (mm) {
        mmgrab(mm);
    }

    // pr_info("Before going into usm, iterating and remove page tables:\n");
    // iterate_hash_table();

    pr_info("User madvise address %ld, size %ld, vma %px, mm %px\n",\
            start, end - start, vma, mm);

	for (cur_addr = start; cur_addr < end; cur_addr += PAGE_SIZE) {
        // down_read(&mm->mmap_sem);
        pr_info("counter: %d\n", ++counter);
        pr_info("checking cur_addr %ld\n", cur_addr);
        /* get the struct page at address cur_addr */
        cur_page = follow_page(vma, cur_addr, FOLL_GET);
        pr_info("cur_page: %px, mapcount %d, refcount %d\n", cur_page, page_mapcount(cur_page), page_ref_count(cur_page));
        // up_read(&mm->mmap_sem);
        if (cur_page == NULL) {
            pr_info("usm: no present page at address %ld\n", cur_addr);
            continue;
        }

        hash_value = get_hash_value(cur_page);
        hash_index = get_hash_index(hash_value);

        /* check if the page is already stored by USM */
        spin_lock(&rmap_hash_lock);
        rmap_bucket = &(rmap_hash_table[cur_addr % nrmap_hash]);
        hlist_for_each_entry_safe(old_rmap_node, n, rmap_bucket, hlist_link) {
            if (old_rmap_node->addr == cur_addr && old_rmap_node->mm == mm) {
                if (old_rmap_node->old_hash_value != hash_value) {
                    hash_del(&old_rmap_node->hlist_link);
                    spin_lock(&page_hash_lock);
                    page_bucket = &(page_hash_table[old_rmap_node->old_hash_value % npage_hash]);
                    hlist_for_each_entry_safe(old_page_node, cur_node, page_bucket, hlist_link) {
                        if (old_page_node->hash_value == old_rmap_node->old_hash_value) {
                            hash_del(&old_page_node->hlist_link);
                            break;
                        }
                    }
                    spin_unlock(&page_hash_lock);
                    spin_unlock(&rmap_hash_lock);
                    goto search_hash_table;
                }
                else {
                    /* the page is already added to USM, no need to do anything */
                    spin_unlock(&rmap_hash_lock);
                    goto next;
                }
            }
        }
        spin_unlock(&rmap_hash_lock);
    search_hash_table:
        search_hash_table(hash_value, hash_index, cur_page, cur_addr, mm, vma, vm_flags);
    next:
        // iterate_hash_table();
        put_page(cur_page);
        pr_info("put_page after follow page for madvised page %px\n", cur_page);
        continue;
    }

    mmdrop(mm);
    return 0;
}

/* search in the mm the pages stored in the hash tables, and delete then */
void usm_exit(struct mm_struct *mm)
{
    struct vm_area_struct *vma = mm->mmap;
    unsigned long addr;
    struct page_node *cur_page_node;
    struct hlist_node *tmp_page_node;
    struct rmap_node *cur_rmap_node;
    struct hlist_node *tmp_rmap_node;
    u64 hash_val;

    if (!test_bit(MMF_VM_SHAREABLE, &mm->flags))
        return;

    pr_info("usm_exit\n");
    while (vma) {
        if (vma->vm_flags & VM_SHAREABLE) {
            addr = vma->vm_start;
            /* iterate over each page in the vma */
            while (addr < vma->vm_end) {
                hlist_for_each_entry_safe(cur_rmap_node, tmp_rmap_node,\
                    &rmap_hash_table[addr % nrmap_hash], hlist_link) {
                    if (cur_rmap_node->mm != mm || cur_rmap_node->addr != addr)
                        continue;
                    hash_val = cur_rmap_node->old_hash_value;
                    hlist_for_each_entry_safe(cur_page_node, tmp_page_node,\
                        &page_hash_table[hash_val % npage_hash], hlist_link) {
                        if (cur_page_node->hash_value != hash_val ||\
                            cur_page_node->mm != mm || cur_page_node->addr != addr)
                            continue;
                        hash_del(&cur_page_node->hlist_link);
                        free_page_node(cur_page_node);
                        break;
                    }
                    hash_del(&cur_rmap_node->hlist_link);
                    free_rmap_node(cur_rmap_node);
                    break;
                }
                addr += PAGE_SIZE;
            }
            vma->vm_flags &= ~VM_SHAREABLE;
        }
        vma = vma->vm_next;
    }
    clear_bit(MMF_VM_SHAREABLE, &mm->flags);
    iterate_hash_table();
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
