#ifndef _ASM_X86_MCE_H
#define _ASM_X86_MCE_H

#include <linux/types.h>
#include <asm/ioctls.h>

/*
 * Machine Check support for x86
 */

/* MCG_CAP register defines */
#define MCG_BANKCNT_MASK	0xff         /* Number of Banks */
#define MCG_CTL_P		(1ULL<<8)    /* MCG_CTL register available */
#define MCG_EXT_P		(1ULL<<9)    /* Extended registers available */
#define MCG_CMCI_P		(1ULL<<10)   /* CMCI supported */
#define MCG_EXT_CNT_MASK	0xff0000     /* Number of Extended registers */
#define MCG_EXT_CNT_SHIFT	16
#define MCG_EXT_CNT(c)		(((c) & MCG_EXT_CNT_MASK) >> MCG_EXT_CNT_SHIFT)
#define MCG_SER_P	 	(1ULL<<24)   /* MCA recovery/new status bits */
#define MCG_ELOG_P		(1ULL<<26)   /* Extended error log supported */

/* MCG_STATUS register defines */
#define MCG_STATUS_RIPV  (1ULL<<0)   /* restart ip valid */
#define MCG_STATUS_EIPV  (1ULL<<1)   /* ip points to correct instruction */
#define MCG_STATUS_MCIP  (1ULL<<2)   /* machine check in progress */

/* MCi_STATUS register defines */
#define MCI_STATUS_VAL   (1ULL<<63)  /* valid error */
#define MCI_STATUS_OVER  (1ULL<<62)  /* previous errors lost */
#define MCI_STATUS_UC    (1ULL<<61)  /* uncorrected error */
#define MCI_STATUS_EN    (1ULL<<60)  /* error enabled */
#define MCI_STATUS_MISCV (1ULL<<59)  /* misc error reg. valid */
#define MCI_STATUS_ADDRV (1ULL<<58)  /* addr reg. valid */
#define MCI_STATUS_PCC   (1ULL<<57)  /* processor context corrupt */
#define MCI_STATUS_S	 (1ULL<<56)  /* Signaled machine check */
#define MCI_STATUS_AR	 (1ULL<<55)  /* Action required */

/*
 * Note that the full MCACOD field of IA32_MCi_STATUS MSR is
 * bits 15:0.  But bit 12 is the 'F' bit, defined for corrected
 * errors to indicate that errors are being filtered by hardware.
 * We should mask out bit 12 when looking for specific signatures
 * of uncorrected errors - so the F bit is deliberately skipped
 * in this #define.
 */
#define MCACOD		  0xefff     /* MCA Error Code */

/* Architecturally defined codes from SDM Vol. 3B Chapter 15 */
#define MCACOD_SCRUB	0x00C0	/* 0xC0-0xCF Memory Scrubbing */
#define MCACOD_SCRUBMSK	0xeff0	/* Skip bit 12 ('F' bit) */
#define MCACOD_L3WB	0x017A	/* L3 Explicit Writeback */
#define MCACOD_DATA	0x0134	/* Data Load */
#define MCACOD_INSTR	0x0150	/* Instruction Fetch */

/* MCi_MISC register defines */
#define MCI_MISC_ADDR_LSB(m)	((m) & 0x3f)
#define MCI_MISC_ADDR_MODE(m)	(((m) >> 6) & 7)
#define  MCI_MISC_ADDR_SEGOFF	0	/* segment offset */
#define  MCI_MISC_ADDR_LINEAR	1	/* linear address */
#define  MCI_MISC_ADDR_PHYS	2	/* physical address */
#define  MCI_MISC_ADDR_MEM	3	/* memory address */
#define  MCI_MISC_ADDR_GENERIC	7	/* generic */

/* CTL2 register defines */
#define MCI_CTL2_CMCI_EN		(1ULL << 30)
#define MCI_CTL2_CMCI_THRESHOLD_MASK	0x7fffULL

#define MCJ_CTX_MASK		3
#define MCJ_CTX(flags)		((flags) & MCJ_CTX_MASK)
#define MCJ_CTX_RANDOM		0    /* inject context: random */
#define MCJ_CTX_PROCESS		1    /* inject context: process */
#define MCJ_CTX_IRQ		2    /* inject context: IRQ */
#define MCJ_NMI_BROADCAST	4    /* do NMI broadcasting */
#define MCJ_EXCEPTION		8    /* raise as exception */

/* Fields are zero when not available */
struct mce {
	__u64 status;
	__u64 misc;
	__u64 addr;
	__u64 mcgstatus;
	__u64 ip;
	__u64 tsc;	/* cpu time stamp counter */
	__u64 time;	/* wall time_t when error was detected */
	__u8  cpuvendor;	/* cpu vendor as encoded in system.h */
	__u8  inject_flags;	/* software inject flags */
	__u16  pad;
	__u32 cpuid;	/* CPUID 1 EAX */
	__u8  cs;		/* code segment */
	__u8  bank;	/* machine check bank */
	__u8  cpu;	/* cpu number; obsolete; use extcpu now */
	__u8  finished;   /* entry is valid */
	__u32 extcpu;	/* linux cpu number that detected the error */
	__u32 socketid;	/* CPU socket ID */
	__u32 apicid;	/* CPU initial apic ID */
	__u64 mcgcap;	/* MCGCAP MSR: machine check capabilities of CPU */
};

/*
 * This structure contains all data related to the MCE log.  Also
 * carries a signature to make it easier to find from external
 * debugging tools.  Each entry is only valid when its finished flag
 * is set.
 */

#define MCE_LOG_LEN 32

struct mce_log {
	char signature[12]; /* "MACHINECHECK" */
	unsigned len;	    /* = MCE_LOG_LEN */
	unsigned next;
	unsigned flags;
	unsigned recordlen;	/* length of struct mce */
	struct mce entry[MCE_LOG_LEN];
};

#define MCE_OVERFLOW 0		/* bit 0 in flags means overflow */

#define MCE_LOG_SIGNATURE	"MACHINECHECK"

#define MCE_GET_RECORD_LEN   _IOR('M', 1, int)
#define MCE_GET_LOG_LEN      _IOR('M', 2, int)
#define MCE_GETCLEAR_FLAGS   _IOR('M', 3, int)

/* Software defined banks */
#define MCE_EXTENDED_BANK	128
#define MCE_THERMAL_BANK	MCE_EXTENDED_BANK + 0

#define K8_MCE_THRESHOLD_BASE      (MCE_EXTENDED_BANK + 1)      /* MCE_AMD */
#define K8_MCE_THRESHOLD_BANK_0    (MCE_THRESHOLD_BASE + 0 * 9)
#define K8_MCE_THRESHOLD_BANK_1    (MCE_THRESHOLD_BASE + 1 * 9)
#define K8_MCE_THRESHOLD_BANK_2    (MCE_THRESHOLD_BASE + 2 * 9)
#define K8_MCE_THRESHOLD_BANK_3    (MCE_THRESHOLD_BASE + 3 * 9)
#define K8_MCE_THRESHOLD_BANK_4    (MCE_THRESHOLD_BASE + 4 * 9)
#define K8_MCE_THRESHOLD_BANK_5    (MCE_THRESHOLD_BASE + 5 * 9)
#define K8_MCE_THRESHOLD_DRAM_ECC  (MCE_THRESHOLD_BANK_4 + 0)


#endif /* _ASM_X86_MCE_H */
