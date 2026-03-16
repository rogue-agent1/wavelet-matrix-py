#!/usr/bin/env python3
"""Wavelet matrix — rank, select, and quantile queries on integer arrays."""

class WaveletMatrix:
    def __init__(self, arr, sigma=None):
        n = len(arr); self.n = n
        self.sigma = sigma or (max(arr)+1 if arr else 1)
        self.lg = max(1, (self.sigma-1).bit_length())
        self.bv = []; self.z = []
        cur = list(arr)
        for d in range(self.lg-1, -1, -1):
            bits = [(v >> d) & 1 for v in cur]
            self.bv.append(bits)
            zeros = [v for v, b in zip(cur, bits) if b == 0]
            ones = [v for v, b in zip(cur, bits) if b == 1]
            self.z.append(len(zeros))
            cur = zeros + ones
        # Build prefix sums for each bit-vector
        self.ps = []
        for bits in self.bv:
            ps = [0]*(n+1)
            for i, b in enumerate(bits): ps[i+1] = ps[i] + b
            self.ps.append(ps)
    def _rank1(self, level, pos):
        return self.ps[level][pos]
    def _rank0(self, level, pos):
        return pos - self._rank1(level, pos)
    def quantile(self, l, r, k):
        """k-th smallest in arr[l:r+1], 0-indexed."""
        for d in range(len(self.bv)):
            zeros_l = self._rank0(d, l); zeros_r = self._rank0(d, r+1)
            cnt0 = zeros_r - zeros_l
            if k < cnt0:
                l = zeros_l; r = zeros_r - 1
            else:
                k -= cnt0
                l = self.z[d] + self._rank1(d, l)
                r = self.z[d] + self._rank1(d, r+1) - 1
        return l  # position in final sorted order... we need original value
    def kth(self, l, r, k):
        """k-th smallest value in arr[l..r], 0-indexed."""
        val = 0
        for d in range(len(self.bv)):
            zeros_l = self._rank0(d, l); zeros_r = self._rank0(d, r+1)
            cnt0 = zeros_r - zeros_l
            bit_pos = self.lg - 1 - d
            if k < cnt0:
                l = zeros_l; r = zeros_r - 1
            else:
                k -= cnt0; val |= (1 << bit_pos)
                l = self.z[d] + self._rank1(d, l)
                r = self.z[d] + self._rank1(d, r+1) - 1
        return val

def main():
    arr = [3,1,4,1,5,9,2,6]
    wm = WaveletMatrix(arr, 10)
    print(f"2nd smallest in [0,7]: {wm.kth(0,7,1)}")
    print(f"Median of [2,5]: {wm.kth(2,5,1)}")

if __name__ == "__main__": main()
