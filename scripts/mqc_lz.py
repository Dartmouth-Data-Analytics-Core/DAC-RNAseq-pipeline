"""Pure-Python LZString decoder (decompressFromBase64).
Replaces the `lzstring` package, which does not install on Python >= 3.12."""

_KEY = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/="
_VAL = {c: i for i, c in enumerate(_KEY)}


def decompress_from_base64(s: str) -> str:
    s = s.strip()
    if not s:
        return ""
    vals = [_VAL[c] for c in s if c in _VAL]
    n = len(vals)
    pos = 0
    val = vals[0]
    reset = 32
    idx = 1

    def read(bits):
        nonlocal val, pos, idx, reset
        res, power, maxp = 0, 1, 1 << bits
        while power != maxp:
            b = val & reset
            reset >>= 1
            if reset == 0:
                reset = 32
                val = vals[idx] if idx < n else 0
                idx += 1
            if b:
                res |= power
            power <<= 1
        return res

    dictionary = {0: 0, 1: 1, 2: 2}
    enlarge_in, dict_size, num_bits = 4, 4, 3
    c = read(2)
    if c == 2:
        return ""
    w = chr(read(8 if c == 0 else 16))
    dictionary[3] = w
    out = [w]
    while True:
        if idx > n + 1:
            return "".join(out)
        c = read(num_bits)
        if c in (0, 1):
            dictionary[dict_size] = chr(read(8 if c == 0 else 16))
            dict_size += 1
            c = dict_size - 1
            enlarge_in -= 1
        elif c == 2:
            return "".join(out)
        if enlarge_in == 0:
            enlarge_in = 1 << num_bits
            num_bits += 1
        if c in dictionary:
            entry = dictionary[c]
        elif c == dict_size:
            entry = w + w[0]
        else:
            raise ValueError("LZString: corrupted data")
        out.append(entry)
        dictionary[dict_size] = w + entry[0]
        dict_size += 1
        enlarge_in -= 1
        w = entry
        if enlarge_in == 0:
            enlarge_in = 1 << num_bits
            num_bits += 1
