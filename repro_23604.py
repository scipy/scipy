"""
Reproduce scipy/scipy#23604: find_peaks_cwt bad SNR estimate.

Shows that the current _filter_ridge_lines uses the CWT value at the
END of the ridge (smallest scale), where values are fading off, instead
of the MAXIMUM value along the ridge as specified in the reference paper.

With the bug: True positives drop off at low min_snr values.
With the fix: True positives persist to much higher min_snr values.
"""
import numpy as np
from scipy.signal import find_peaks_cwt


def lorentz(x, height, width, loc):
    return height / (1 + np.power((x - loc) / width, 2))


rng = np.random.default_rng(seed=50)
x = np.arange(500)
noise = rng.standard_normal(500) * 1e-2
signal = (
    lorentz(x, 0.1, 10, 100)
    + lorentz(x, 1, 15, 250)
    + lorentz(x, 0.5, 5, 400)
    + noise
)

cwt_scales = np.geomspace(1, 20, 16)

# Test with min_snr=1 (should detect peaks easily)
peaks_low = find_peaks_cwt(signal, cwt_scales, min_snr=1, noise_perc=90)
# Test with min_snr=50 (higher SNR threshold)
peaks_high = find_peaks_cwt(signal, cwt_scales, min_snr=50, noise_perc=90)
# Test with min_snr=100
peaks_vhigh = find_peaks_cwt(signal, cwt_scales, min_snr=100, noise_perc=90)

true_peaks = {100, 250, 400}

print("=== Reproducing scipy/scipy#23604 ===")
print(f"True peaks at indices: {true_peaks}")
print()
print(f"min_snr=1:  found {peaks_low}")
tp_low = set(peaks_low) & true_peaks
fp_low = set(peaks_low) - true_peaks
print(f"  True positives: {tp_low}")
print(f"  False positives: {fp_low}")
print()
print(f"min_snr=50: found {peaks_high}")
tp_high = set(peaks_high) & true_peaks
fp_high = set(peaks_high) - true_peaks
print(f"  True positives: {tp_high}")
print(f"  False positives: {fp_high}")
print()
print(f"min_snr=100: found {peaks_vhigh}")
tp_vhigh = set(peaks_vhigh) & true_peaks
fp_vhigh = set(peaks_vhigh) - true_peaks
print(f"  True positives: {tp_vhigh}")
print(f"  False positives: {fp_vhigh}")
print()

# Find the min_snr at which we lose the first true peak
print("=== Finding the SNR threshold where true peaks are lost ===")
for snrval in np.linspace(1, 200, 200):
    peaks = find_peaks_cwt(signal, cwt_scales, min_snr=snrval, noise_perc=90)
    tps = set(peaks) & true_peaks
    if len(tps) < 3:
        print(f"At min_snr={snrval:.1f}: lost {true_peaks - tps}, keeping {tps}")
        break
else:
    print("All true peaks detected up to min_snr=200")
