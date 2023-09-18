import functools
import re
from typing import Dict, List, Sequence

DECODE_DEFAULT_CHARS = ";/?:@&=+$,#"
DECODE_COMPONENT_CHARS = ""

decode_cache: Dict[str, List[str]] = {}


def get_decode_cache(exclude: str) -> Sequence[str]:
    if exclude in decode_cache:
        return decode_cache[exclude]

    cache: List[str] = []
    decode_cache[exclude] = cache

    for i in range(128):
        ch = chr(i)
        cache.append(ch)

    for i in range(len(exclude)):
        ch_code = ord(exclude[i])
        cache[ch_code] = "%" + ("0" + hex(ch_code)[2:].upper())[-2:]

    return cache


# Decode percent-encoded string.
#
def decode(string: str, exclude: str = DECODE_DEFAULT_CHARS) -> str:
    cache = get_decode_cache(exclude)
    repl_func = functools.partial(repl_func_with_cache, cache=cache)
    return re.sub(r"(%[a-f0-9]{2})+", repl_func, string, flags=re.IGNORECASE)


def repl_func_with_cache(match: "re.Match", cache: Sequence[str]) -> str:
    seq = match.group()
    result = ""

    i = 0
    l = len(seq)  # noqa: E741
    while i < l:
        b1 = int(seq[i + 1 : i + 3], 16)

        if b1 < 0x80:
            result += cache[b1]
            i += 3  # emulate JS for loop statement3
            continue

        if (b1 & 0xE0) == 0xC0 and (i + 3 < l):
            # 110xxxxx 10xxxxxx
            b2 = int(seq[i + 4 : i + 6], 16)

            if (b2 & 0xC0) == 0x80:
                chr_ = ((b1 << 6) & 0x7C0) | (b2 & 0x3F)

                if chr_ < 0x80:
                    result += "\ufffd\ufffd"
                else:
                    result += chr(chr_)

                i += 3
                i += 3  # emulate JS for loop statement3
                continue

        if (b1 & 0xF0) == 0xE0 and (i + 6 < l):
            # 1110xxxx 10xxxxxx 10xxxxxx
            b2 = int(seq[i + 4 : i + 6], 16)
            b3 = int(seq[i + 7 : i + 9], 16)

            if (b2 & 0xC0) == 0x80 and (b3 & 0xC0) == 0x80:
                chr_ = ((b1 << 12) & 0xF000) | ((b2 << 6) & 0xFC0) | (b3 & 0x3F)

                if chr_ < 0x800 or (chr_ >= 0xD800 and chr_ <= 0xDFFF):
                    result += "\ufffd\ufffd\ufffd"
                else:
                    result += chr(chr_)

                i += 6
                i += 3  # emulate JS for loop statement3
                continue

        if (b1 & 0xF8) == 0xF0 and (i + 9 < l):
            # 111110xx 10xxxxxx 10xxxxxx 10xxxxxx
            b2 = int(seq[i + 4 : i + 6], 16)
            b3 = int(seq[i + 7 : i + 9], 16)
            b4 = int(seq[i + 10 : i + 12], 16)

            if (b2 & 0xC0) == 0x80 and (b3 & 0xC0) == 0x80 and (b4 & 0xC0) == 0x80:
                chr_ = (
                    ((b1 << 18) & 0x1C0000)
                    | ((b2 << 12) & 0x3F000)
                    | ((b3 << 6) & 0xFC0)
                    | (b4 & 0x3F)
                )

                if chr_ < 0x10000 or chr_ > 0x10FFFF:
                    result += "\ufffd\ufffd\ufffd\ufffd"
                else:
                    chr_ -= 0x10000
                    result += chr(0xD800 + (chr_ >> 10)) + chr(0xDC00 + (chr_ & 0x3FF))

                i += 9
                i += 3  # emulate JS for loop statement3
                continue

        result += "\ufffd"
        i += 3  # emulate JS for loop statement3

    return result
