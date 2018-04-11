"""
Functions that load datasets.
"""

from __future__ import division, print_function, absolute_import

from numpy import array, frombuffer, load


def ascent():
    """
    Get an 8-bit grayscale 512 x 512 image for easy use in demos.

    The image is derived from accent-to-the-top.jpg at
    http://www.public-domain-image.com/people-public-domain-images-pictures/

    Parameters
    ----------
    None

    Returns
    -------
    ascent : ndarray
       convenient image to use for testing and demonstration

    Examples
    --------
    >>> import scipy.datasets
    >>> ascent = scipy.datasets.ascent()
    >>> ascent.shape
    (512, 512)
    >>> ascent.max()
    255

    >>> import matplotlib.pyplot as plt
    >>> plt.gray()
    >>> plt.imshow(ascent)
    >>> plt.show()

    """
    import pickle
    import os
    fname = os.path.join(os.path.dirname(__file__), '_data', 'ascent.dat')
    with open(fname, 'rb') as f:
        ascent = array(pickle.load(f))
    return ascent


def face(gray=False):
    """
    Get a 1024 x 768, color image of a raccoon face.

    raccoon-procyon-lotor.jpg at http://www.public-domain-image.com

    Parameters
    ----------
    gray : bool, optional
        If True return 8-bit grey-scale image, otherwise return a color image

    Returns
    -------
    face : ndarray
        image of a racoon face

    Examples
    --------
    >>> import scipy.datasets
    >>> face = scipy.datasets.face()
    >>> face.shape
    (768, 1024, 3)
    >>> face.max()
    255
    >>> face.dtype
    dtype('uint8')

    >>> import matplotlib.pyplot as plt
    >>> plt.gray()
    >>> plt.imshow(face)
    >>> plt.show()

    """
    import bz2
    import os
    with open(os.path.join(os.path.dirname(__file__), '_data', 'face.dat'),
              'rb') as f:
        rawdata = f.read()
    data = bz2.decompress(rawdata)
    face = frombuffer(data, dtype='uint8')
    face.shape = (768, 1024, 3)
    if gray is True:
        face = (0.21 * face[:, :, 0] + 0.71 * face[:, :, 1] +
                0.07 * face[:, :, 2]).astype('uint8')
    return face


def electrocardiogram():
    """
    Load an electrocardiogram as an example for a one-dimensional signal.

    The returned signal is a 5 minute long electrocardiogram (ECG), a medical
    recording of the heart's electrical activity, sampled at 360 Hz.

    Returns
    -------
    ecg : ndarray
        The electrocardiogram in millivolt (mV) sampled at 360 Hz.

    Notes
    -----
    The provided signal is an excerpt (19:35 to 24:35) from the `record 208`_
    (lead MLII) provided by the MIT-BIH Arrhythmia Database [1]_ on
    PhysioNet [2]_. The excerpt includes noise induced artifacts, typical
    heartbeats as well as pathological changes.

    .. _record 208: https://physionet.org/physiobank/database/html/mitdbdir/records.htm#208

    .. versionadded:: 1.1.0

    References
    ----------
    .. [1] Moody GB, Mark RG. The impact of the MIT-BIH Arrhythmia Database.
           IEEE Eng in Med and Biol 20(3):45-50 (May-June 2001).
           (PMID: 11446209); https://doi.org/10.13026/C2F305
    .. [2] Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PCh,
           Mark RG, Mietus JE, Moody GB, Peng C-K, Stanley HE. PhysioBank,
           PhysioToolkit, and PhysioNet: Components of a New Research Resource
           for Complex Physiologic Signals. Circulation 101(23):e215-e220;
           https://doi.org/10.1161/01.CIR.101.23.e215

    Examples
    --------
    >>> from scipy.datasets import electrocardiogram
    >>> ecg = electrocardiogram()
    >>> ecg
    array([-0.245, -0.215, -0.185, ..., -0.405, -0.395, -0.385])
    >>> ecg.shape, ecg.mean(), ecg.std()
    ((108000,), -0.16510875, 0.5992473991177294)

    As stated the signal features several areas with a different morphology.
    E.g. the first few seconds show the electrical activity of a heart in
    normal sinus rhythm as seen below.

    >>> import matplotlib.pyplot as plt
    >>> fs = 360
    >>> time = np.arange(ecg.size) / fs
    >>> plt.plot(time, ecg)
    >>> plt.xlabel("time in s")
    >>> plt.ylabel("ECG in mV")
    >>> plt.xlim(9, 10.2)
    >>> plt.ylim(-1, 1.5)
    >>> plt.show()

    After second 16 however, the first premature ventricular contractions, also
    called extrasystoles, appear. These have a different morphology compared to
    typical heartbeats. The difference can easily be observed in the following
    plot.

    >>> plt.plot(time, ecg)
    >>> plt.xlabel("time in s")
    >>> plt.ylabel("ECG in mV")
    >>> plt.xlim(46.5, 50)
    >>> plt.ylim(-2, 1.5)
    >>> plt.show()

    At several points large artifacts disturb the recording, e.g.:

    >>> plt.plot(time, ecg)
    >>> plt.xlabel("time in s")
    >>> plt.ylabel("ECG in mV")
    >>> plt.xlim(207, 215)
    >>> plt.ylim(-2, 3.5)
    >>> plt.show()

    Finally, examining the power spectrum reveals that most of the biosignal is
    made up of lower frequencies. At 60 Hz the noise induced by the mains
    electricity can be clearly observed.

    >>> from scipy.signal import welch
    >>> f, Pxx = welch(ecg, fs=fs, nperseg=2048, scaling="spectrum")
    >>> plt.semilogy(f, Pxx)
    >>> plt.xlabel("Frequency in Hz")
    >>> plt.ylabel("Power spectrum of the ECG in mV**2")
    >>> plt.xlim(f[[0, -1]])
    >>> plt.show()
    """
    import os
    file_path = os.path.join(os.path.dirname(__file__), '_data', 'ecg.dat')
    with load(file_path) as file:
        ecg = file["ecg"].astype(int)  # np.uint16 -> int
    # Convert raw output of ADC to mV: (ecg - adc_zero) / adc_gain
    ecg = (ecg - 1024) / 200.0
    return ecg
