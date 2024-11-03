import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

from skbio.stats.distance import permanova
from skbio.stats.ordination import pcoa
from skbio.diversity import beta_diversity

from Model.get_lineage import *
from Model.modificator import otu_table, tax_table

from PyQt6.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QHBoxLayout
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=3, height=2, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)

class MicrobiomeDataAnalyzer:

    def __init__(self, OTU_table: str = "", Taxa_table: str = "", Metadata: str = "") -> None:
        self.OTU_table = OTU_table
        self.Taxa_table = Taxa_table
        self.Metadata = Metadata

    def beta_diversity(self):
        dropped_df = self.OTU_table.dropna()
        otumat_transposed = dropped_df.T
        data = otumat_transposed.values
        ids = otumat_transposed.index
        b_div = beta_diversity(metric="braycurtis", counts=data, ids=ids)
        return b_div

    def plot_rank(self, rank, canvas):
        taxa_mapping = self.Taxa_table[rank]
        rank_df = self.OTU_table.groupby(taxa_mapping, axis=0).sum().T

        rank_df.plot(kind='bar', stacked=True, ax=canvas.axes)

        canvas.axes.set_ylabel('Abundance', fontsize=20)
        canvas.axes.tick_params(axis='x', labelsize=14, rotation=45)
        canvas.axes.tick_params(axis='y', labelsize=14)
        canvas.axes.legend(title=rank, fontsize='10', title_fontsize='12', loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)

        # Adjust layout to fit everything nicely
        canvas.figure.tight_layout()

    def plot_top(self, top, canvas):
        total_counts = self.OTU_table.sum(axis=1)
        top_otus = total_counts.nlargest(top).index

        top_otus_tax_rank = {}
        for otu in top_otus:
            tax_info = self.Taxa_table.loc[otu]
            specific_rank = tax_info.dropna().iloc[-1] if not tax_info.dropna().empty else 'Unknown'
            top_otus_tax_rank[otu] = specific_rank

        filtered_df = self.OTU_table.loc[top_otus]
        df_t = filtered_df.T

        df_t.plot(kind='bar', stacked=True, ax=canvas.axes)

        canvas.axes.set_title(f'Top {top} taxa by read count', fontsize=22)
        canvas.axes.set_ylabel('Read Count', fontsize=20)
        canvas.axes.tick_params(axis='x', labelsize=14, rotation=45)
        canvas.axes.tick_params(axis='y', labelsize=14)

        legend_labels = [top_otus_tax_rank[otu] for otu in top_otus]
        canvas.axes.legend(legend_labels, fontsize='10', title_fontsize='12', loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)

        # Adjust layout to fit everything nicely
        canvas.figure.tight_layout()

    def plot_pcoa(self, canvas):
        distance = self.beta_diversity()
        pcoa_results = pcoa(distance)
        pcoa_samples = pcoa_results.samples

        groups = np.array(self.Metadata['Group'])
        treatment_day = np.array(self.Metadata['Property'])
        pcoa_samples.insert(0, "Group", groups)
        pcoa_samples.insert(1, "Property", treatment_day)

        unique_groups = pcoa_samples['Group'].unique()
        colors = plt.get_cmap('tab10')(range(len(unique_groups)))
        color_map = dict(zip(unique_groups, colors))

        seen_groups = set()

        for idx, row in pcoa_samples.iterrows():
            group = row['Group']
            pc1 = row['PC1']
            pc2 = row['PC2']
            property = row['Property']

            if group not in seen_groups:
                canvas.axes.scatter(pc1, pc2, color=color_map[group], label=group, alpha=0.7)
                seen_groups.add(group)
            else:
                canvas.axes.scatter(pc1, pc2, color=color_map[group], alpha=0.7)

            canvas.axes.text(pc1, pc2, property, fontsize=8)

        canvas.axes.set_xlabel('PC1', fontsize=14)
        canvas.axes.set_ylabel('PC2', fontsize=14)
        canvas.axes.set_title('PCoA of Beta Diversity', fontsize=16)
        canvas.axes.legend(loc='upper left', bbox_to_anchor=(1, 1))

class MainWindow(QMainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        self.setWindowTitle("Microbiome Data Analyzer")
        self.setGeometry(50, 50, 1200, 800)

        layout = QHBoxLayout()

        # Create the maptlotlib FigureCanvas object,
        # which defines a single set of axes as self.axes.
        sc_rank = MplCanvas(self, width=3, height=2, dpi=100)
        sc_top = MplCanvas(self, width=3, height=2, dpi=100)
        sc_pcoa = MplCanvas(self, width=3, height=2, dpi=100)

        layout.addWidget(sc_rank)
        layout.addWidget(sc_top)
        layout.addWidget(sc_pcoa)


        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

        # Example usage
        otumat_alice_bob = otu_table(r'C:\Users\egoro\Desktop\MVC\Model\alice_bob.csv')  # Replace with actual OTU table data
        taxmat_alice_bob = tax_table(r'C:\Users\egoro\Desktop\MVC\Model\alice_bob.csv')  # Replace with actual Taxa table data
        metadata_dict_alice_bob = {
            'SampleID': ['Alice00-1mio.daa', 'Alice01-1mio.daa',
                         'Alice03-1mio.daa', 'Alice06-1mio.daa',
                         'Alice08-1mio.daa', 'Alice34-1mio.daa',
                         'Bob00-1mio.daa', 'Bob01-1mio.daa',
                         'Bob03-1mio.daa', 'Bob06-1mio.daa',
                         'Bob08-1mio.daa', 'Bob34-1mio.daa'],
            'Group': ['Alice', 'Alice', 'Alice',
                      'Alice', 'Alice', 'Alice',
                      'Bob', 'Bob', 'Bob',
                      'Bob', 'Bob', 'Bob'],
            'Property': ['0-', '1+', '3+', '6+', '8-', '34-',
                         '0-', '1+', '3+', '6+', '8-', '34-']
        }
        metadata = pd.DataFrame(metadata_dict_alice_bob)
        analyzer = MicrobiomeDataAnalyzer(otumat_alice_bob, taxmat_alice_bob, metadata)

        analyzer.plot_top(5, sc_top)
        analyzer.plot_rank("Genus", sc_rank)
        analyzer.plot_pcoa(sc_pcoa)

if __name__ == "__main__":
    app = QApplication([])
    window = MainWindow()
    window.show()
    app.exec()