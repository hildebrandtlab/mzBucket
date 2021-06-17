import numpy as np
import sys
from itertools import groupby
from functools import reduce
from tqdm import tqdm


def vectorized_window(window, window_length=4, resolution=2):
    """
    :param: window: a window to be vectorized
    :param: windowLength: the length windows were split into
    :param: resolution: number of decimals kept in binned spectrum
    :return: a vectorized window as numpy array
    """
    factor = np.power(10, resolution)
    vectorized = [(int((mz + (window_length / 2.0) * window[1]) * factor)
                   - factor * window[0] * window_length, i) for (mz, i) in window[2]]

    binned_vec_summed = []

    for key, values in groupby(vectorized, lambda x: x[0]):
        binned_vec_summed.append((key, reduce(lambda x, y: x + y, [v[1] for v in list(values)])))

    return window[0], window[1], binned_vec_summed


class Spectrum:

    def __init__(self, mz, intensity, rt_index=0, dt_index=0):
        self.mz = mz
        self.intensity = intensity
        self.rt_index = rt_index
        self.dt_index = dt_index

    def bin_to_resolution(self, resolution=2, min_intensity=150):
        """
        bins a mass spectrum to a finite resolution by rounding to a given precision
        and summing intensities by resulting mz bins
        :param: resolution: number of decimals to keep
        :param: min_intensity: a minimum intensity a peak should have to be considered interesting
        :return: a binned vector as list of tuples [(binned_mz, i)]
        """
        binned_vec = [(np.round(mz, resolution), i) for (mz, i) in list(zip(self.mz, self.intensity))]
        binned_vec = list(filter(lambda x: x[1] > min_intensity, binned_vec))

        binned_vec_summed = []

        for key, values in groupby(binned_vec, lambda x: x[0]):
            binned_vec_summed.append((key, reduce(lambda x, y: x + y, [v[1] for v in list(values)])))

        return binned_vec_summed

    def create_overlapping_windows(self, resolution=2, window_length=4, min_intensity=150):
        """
        splits a mass spectrum with binned resolution into overlapping windows
        :param: mz_i_binned: an binned spectrum of the form [(mz_binned, i)]
        :param: resolution: number of decimals kept in binned spectrum
        :param: window_length: the length windows should be split into
        """
        even_list, odd_list = [], []

        mz_i_binned = self.bin_to_resolution(resolution=resolution, min_intensity=min_intensity)

        # calculate binning keys without offset...
        even = [(int(np.floor(mz / window_length)), mz, i) for (mz, i) in mz_i_binned]

        # ... and with offset
        odd = [(int(np.floor((mz + (window_length / 2.0)) / window_length)), mz, i) for (mz, i) in mz_i_binned]

        # now the grouping part, grouping even...
        for key, values in groupby(even, lambda x: x[0]):
            even_list.append((key, 0, [(mz, i) for (key, mz, i) in values]))

        # ...and odd
        for key, values in groupby(odd, lambda x: x[0]):
            odd_list.append((key, 1, [(mz, i) for (key, mz, i) in values]))

        # concat and return
        return even_list + odd_list

    def sparse_windows_to_dense(self, window_length=4, resolution=2, min_intensity=150):
        """
        :param: windowLength: the length windows were split into
        :param: resolution: number of decimals kept in binned spectrum
        :param: min_intensity: a minimum intensity a peak needs to pe considered interesting for prediction
        :return: a dense numpy array of a given sparse vectorized window
        """
        binned_windows = self.create_overlapping_windows(resolution, window_length, min_intensity)

        all_out = []

        for window in tqdm(binned_windows, desc='creating vectors'):

            if np.max([i for (mz, i) in window[2]]) > 0:
                vec_window = vectorized_window(window, window_length, resolution)
                mz_bin, oe = vec_window[0], vec_window[1]

                col_indices = [mz for (mz, i) in vec_window[2]]
                channel_indices = np.repeat(0, len(col_indices))
                values = [i for (mz, i) in vec_window[2]]

                indices = np.array(list(zip(col_indices, channel_indices)))
                dense = tf.sparse.to_dense(tf.sparse.SparseTensor(indices, values, (1000, 1))).numpy()
                all_out.append((mz_bin, oe, dense))

        return all_out

    def search_isotopic_patterns(self, path_to_model, resolution=2, window_length=4, min_intensity=150):
        """
        :param: path_to_model: path to the tensorflow model
        :param: resolution: number of decimals kept in binned spectrum
        :param: windowLength: the length windows were split into
        :param: min_intensity: a minimum intensity a peak needs to pe considered interesting for prediction
        :return: a set of predicted windows
        """

        # check given parameters
        if window_length > 10.0:
            print("Current version of PeakNet only supports window_lengths up to 10 Da, aborting...")
            sys.exit(1)
        if resolution != 2:
            print("Current version of PeakNet only supports a resolution of 10^-2 Da, aborting....")

        # create dense window collection
        dense_windows = self.sparse_windows_to_dense(window_length, resolution, min_intensity)

        # extract vectors from that collection and normalize for prediction
        window_vec = np.array([w[2] for w in dense_windows])
        window_vec = np.squeeze(window_vec) / np.sum(np.squeeze(window_vec), axis=1)[:, None]

        # instantiate the deep neural network for prediction
        PEAKNET = tf.keras.models.load_model(path_to_model, custom_objects={'loss_charge': 0, 'loss_label': 0})

        # predict charge state and label of all windows that are part of the spectrum
        charges, labels = PEAKNET.predict(np.expand_dims(window_vec, 2))

        # create more readable representation of the outputs
        labels = np.squeeze(np.round(labels))
        charges = np.argmax(charges, axis=1) + 1

        # zip labels and charges onto all windows of the given spectrum
        with_labels = []
        for i, window in enumerate(dense_windows):
            with_labels.append((labels[i], charges[i], window[0], window[1], window[2]))

        # filter to just keep windows that were predicted to contain an iso-pattern and return
        return list(filter(lambda x: x[0] == 1, with_labels))
