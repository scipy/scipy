#!/usr/bin/python
"""
Detect which x86 CPU extension instructions the given scipy install uses.
This file can be used in the release process to check that the nosse installer
does not contain SSE instructions. This has happened before, see for example
ticket #1170.

Is meant to be run on OS X with Wine. Make sure objdump.exe is installed.

See also tools/win32build/misc/x86analysis.py in numpy for a similar script
that checks a single file.

"""
import subprocess
import sys
import os
from optparse import OptionParser


OBJDUMP = os.environ['HOME'] + '/.wine/drive_c/MinGW/bin/objdump.exe'
SCIPY_PY25 = os.environ['HOME'] + '/.wine/drive_c/Python25/Lib/site-packages/scipy/'
SCIPY_PY26 = os.environ['HOME'] + '/.wine/drive_c/Python26/Lib/site-packages/scipy/'
SCIPY_PY27 = os.environ['HOME'] + '/.wine/drive_c/Python27/Lib/site-packages/scipy/'
SCIPY_PY31 = os.environ['HOME'] + '/.wine/drive_c/Python31/Lib/site-packages/scipy/'
NUMPY_PY25 = os.environ['HOME'] + '/.wine/drive_c/Python25/Lib/site-packages/numpy/'
NUMPY_PY26 = os.environ['HOME'] + '/.wine/drive_c/Python26/Lib/site-packages/numpy/'
NUMPY_PY27 = os.environ['HOME'] + '/.wine/drive_c/Python27/Lib/site-packages/numpy/'
NUMPY_PY31 = os.environ['HOME'] + '/.wine/drive_c/Python31/Lib/site-packages/numpy/'
SSE3_LIBS = os.environ['HOME'] + '/.wine/drive_c/local/lib/yop/sse3'
SSE2_LIBS = os.environ['HOME'] + '/.wine/drive_c/local/lib/yop/sse2'
NOSSE_LIBS = os.environ['HOME'] + '/.wine/drive_c/local/lib/yop/nosse'

# The install to check
basepath = SCIPY_PY25


def main():
    # a set of all unique CPU extension codes found
    allcodes = set()
    # walk the SciPy tree and check all binary files
    for root, dirs, files in os.walk(basepath):
        for fl in files:
            if os.path.splitext(fl)[1] in ['.a', '.pyd', '.so']:
                full_fpath = os.path.join(root, fl)
                codes = single_file_checkext(full_fpath)
                for code in codes:
                    allcodes.add(code)
    write_summary(allcodes)


def single_file_checkext(fname, striproot=True):
    if striproot:
        sys.stdout.write('%s: ' % fname.replace(basepath, ''))
    else:
        sys.stdout.write('%s: ' % fname)
    sys.stdout.flush()
    codes = process(path_as_windows(fname))
    sys.stdout.write(" ".join(codes))
    sys.stdout.write("\n")
    return codes


def path_as_windows(fpath):
    """Return the file path as Wine expects."""
    winepath = 'C:\\' + fpath.split('drive_c')[1]
    return winepath


def write_summary(allcodes):
    """Write a summary of all found codes to stdout."""
    print """\n
----------------------------------------------------------------------------
Checked all binary files for CPU extension codes. Found the following codes:"""
    for code in allcodes:
        print code
    print """
----------------------------------------------------------------------------
"""


def process(fn):
    p = subprocess.Popen(['wine', OBJDUMP, '-d', fn], stdout=subprocess.PIPE)
    codes = {}
    for line in p.stdout:
        r = line.split("\t")
        if len(r) != 3:
            continue
        instr = r[2].split()[0].lower()
        if instr in INSTRS:
            codes[INSTRS[instr]] = True
            print instr
    codes = codes.keys()
    codes.sort()
    return codes

#------------------------------------------------------------------------------
# Instruction lists
#------------------------------------------------------------------------------

# x86
EXTS_x86 = dict(
    _486='bswap cmpxch cpuid invd invlpg wbinvd xadd',
    pentium='cmpxchg8b rdmsr rdtsc wrmsr',
    pentium_mmx='rdpmc',
    pentium_pro='cmova cmovae cmovb cmovbe cmovc cmove cmovg cmovge cmovl cmovle cmovna cmovnae cmovnb cmovnbe cmovnc cmovne cmovng cmovnge cmovnl cmovnle cmovno cmovnp cmovns cmovnz cmovo cmovp cmovpe cmovpo cmovs cmovz sysenter sysexit rdpmc ud2',
    amd_k6_2='syscall sysret',
    sse='maskmovq movntps movntq prefetch0 prefetch1 prefetch2 prefetchnta sfence',
    sse2='clflush lfence maskmovdqu mfence movntdq movnti movntpd pause',
    sse3='lddqu',
    sse3_intel='monitor mwait',
    intel_vt='vmptrld vmptrst vmclear vmread vmwrite vmcall vmlaunch vmresume vmxoff vmxon',
    amd_v='clgi skinit stgi vmload vmmcall vmrun vmsave',
    x86_64='cmpxchg16b rdtscp',
    sse4a='lzcnt popcnt',
)


# x87
EXTS_x87 = dict(
    pentium_pro='fcmovb, fcmovbe, fcmove, fcmovnb, fcmovnbe, fcmovne, fcmovnu, fcmovu fcomi fcomip fucomi fucomip',
    sse='fxrstor fxsave',
    sse3='fisttp',
    undocumented='ffreep',
)

# SIMD
EXTS_simd = dict(
    mmx='emms movd movq packssdw packsswb packuswb paddb paddd paddsb paddsw paddusb paddusw paddw pand pandn pcmpeqb pcmpeqd pcmpeqw pcmpgtb pcmpgtd pcmpgtw pmaddwd pmulhw pmullw por pslld psllq psllw psrad psraw psrld psrlq psrlw psubb psubd psubsb psubsw psubusb psubusw psubw punpckhbw punpckhdq punpckhwd punpcklbw punpckldq punpcklwd pxor',
    emmx='paveb paddsiw pmagw pdistib psubsiw pmvzb pmulhrw pmvnzb pmvlzb pmvgezb pmulhriw pmachriw',
    _3dnow='femms pavgusb pf2id pfacc pfadd pfcmpeq pfcmpge pfcmpgt pfmax pfmin pfmul pfrcp pfrcpit1 pfrcpit2 pfrsqit1 pfrsqrt pfsub pfsubr pi2fd pmulhrw prefetch prefetchw',
    _3dnowplus='pf2iw pfnacc pfpnacc pi2fw pswapd',
    _3dnowplus_geodegx='pfrsqrtv pfrcpv',
    sse='addps addss cmpps cmpss comiss cvtpi2ps cvtps2pi cvtsi2ss cvtss2si cvttps2pi cvttss2si divps divss ldmxcsr maxps maxss minps minss movaps movhlps movhps movlhps movlps movmskps movntps movss movups mulps mulss rcpps rcpss rsqrtps rsqrtss shufps sqrtps sqrtss stmxcsr subps subss ucomiss unpckhps unpcklps andnps andps orps pavgb pavgw pextrw pinsrw pmaxsw pmaxub pminsw pminub pmovmskb pmulhuw psadbw pshufw xorps',
    sse2='addpd addsd andnpd andpd cmppd cmpsd comisd cvtdq2pd cvtdq2ps cvtpd2dq cvtpd2pi cvtpd2ps cvtpi2pd cvtps2dq cvtps2pd cvtsd2si cvtsd2ss cvtsi2sd cvtss2sd cvttpd2dq cvttpd2pi cvtps2dq cvttsd2si divpd divsd maxpd maxsd minpd minsd movapd movhpd movlpd movmskpd movsd movupd mulpd mulsd orpd shufpd sqrtpd sqrtsd subpd subsd ucomisd unpckhpd unpcklpd xorpd movdq2q movdqa movdqu movq2dq paddq psubq pmuludq pshufhw pshuflw pshufd pslldq psrldq punpckhqdq punpcklqdq',
    sse3='addsubpd addsubps haddpd haddps hsubpd hsubps movddup movshdup movsldup',
    ssse3='psignw psignd psignb pshufb pmulhrsw pmaddubsw phsubw phsubsw phsubd phaddw phaddsw phaddd palignr pabsw pabsd pabsb',
    sse4_1='mpsadbw phminposuw pmulld pmuldq     dpps dppd     blendps blendpd blendvps blendvpd pblendvb pblendw     pminsb pmaxsb pminuw pmaxuw pminud pmaxud pminsd pmaxsd     roundps roundss roundpd roundsd     insertps pinsrb pinsrd/pinsrq extractps pextrb pextrw pextrd/pextrq     pmovsxbw pmovzxbw pmovsxbd pmovzxbd pmovsxbq pmovzxbq pmovsxwd pmovzxwd pmovsxwq pmovzxwq pmovsxdq pmovzxdq     ptest     pcmpeqq     packusdw     movntdqa',
    sse4a='extrq insertq movntsd movntss',
    sse4_2='crc32     pcmpestri     pcmpestrm     pcmpistri     pcmpistrm     pcmpgtq',
    fma='vfmaddpd vfmaddps vfmaddsd vfmaddss vfmaddsubpd vfmaddsubps vfmsubaddpd vfmsubaddps vfmsubpd vfmsubps vfmsubsd vfmsubss vfnmaddpd vfnmaddps vfnmaddsd vfnmadss vfnmsubpd vfnmsubps vfnmsubsd vfnmsubss',
)

INSTRS = dict()
for ext in [EXTS_x86, EXTS_x87, EXTS_simd]:
    for key, value in ext.items():
        if key.startswith('_'):
            key = key[1:]
        for v in value.split():
            INSTRS[v] = key

#------------------------------------------------------------------------------

if __name__ == "__main__": main()
