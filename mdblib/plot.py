import matplotlib.pyplot as plt
import matplotlib.patches as patches
import decimal
from collections.abc import Sequence
import numbers
from collections import OrderedDict

from .consensus import feature_desc


class linspace(Sequence):
    """linspace(start, stop, num) -> linspace object

    Return a virtual sequence of num numbers from start to stop (inclusive).

    If you need a half-open range, use linspace(start, stop, num+1)[:-1].
    """

    def __init__(self, start, stop, num):
        if not isinstance(num, numbers.Integral) or num <= 1:
            raise ValueError('num must be an integer > 1')
        self.start, self.stop, self.num = start, stop, num
        self.step = (stop - start) / (num - 1)

    def __len__(self):
        return self.num

    def __getitem__(self, i):
        if isinstance(i, slice):
            return [self[x] for x in range(*i.indices(len(self)))]
        if i < 0:
            i = self.num + i
        if i >= self.num:
            raise IndexError('linspace object index out of range')
        if i == self.num - 1:
            return self.stop
        return self.start + i * self.step

    def __repr__(self):
        return '{}({}, {}, {})'.format(type(self).__name__,
                                       self.start, self.stop, self.num)

    def __eq__(self, other):
        if not isinstance(other, linspace):
            return False
        return ((self.start, self.stop, self.num) == (other.start, other.stop, other.num))

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash((type(self), self.start, self.stop, self.num))


def drange(x, y, jump):
    x = decimal.Decimal(x)
    while x < y:
        yield float(x)
        x += decimal.Decimal(jump)


class Plot():
    def __init__(self, consensus, seq):
        feature_order = {y: x for x, y in feature_desc.items()}
        cmap = plt.get_cmap('plasma')
        regs = consensus.prediction.regions

        fig, axes = plt.subplots(2, 1, figsize=(10, 3), sharex='col',
                                 gridspec_kw={'height_ratios': [1, 2]})

        ax = axes[0]
        ax.set_xlim(0, len(seq))
        ax.set_ylim(0, 1.25)

        regbounds = list()
        for reg in regs:
            regbounds.extend(reg[:-1])
            rect = patches.Rectangle((reg[0], 0.75), reg[1] - reg[0] + 1, 0.25)
            ax.add_patch(rect)

        colors = cmap(linspace(0, 0.8, len(feature_desc)))

        for reg in consensus.enriched_regions:
            fdesc = reg[2]
            fcolor_index = feature_order.get(fdesc)
            if fcolor_index:
                fcolor_index = int(fcolor_index) - 1
                rect = patches.Rectangle((reg[0], 0.25), reg[1] - reg[0] + 1, 0.25,
                                         facecolor=colors[fcolor_index], edgecolor='w',
                                         label=fdesc)
                ax.add_patch(rect)

        handles, labels = ax.get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc='upper left', bbox_to_anchor=(1, 1.07))

        ax = axes[1]
        ax.plot(consensus.prediction.scores, linewidth=0.7)
        # ax.fill_between(range(1, len(consensus.prediction.scores) + 1), consensus.prediction.scores, where=[True if y >= 5/8 else False for y in consensus.prediction.scores])
        for coord in regbounds:
            ax.axvline(coord - 1, ls='--', c='k', alpha=0.5, linewidth=0.5)
        ax.axhline(5/8, color='red', linewidth=0.5)
        plt.tight_layout()
        plt.show()